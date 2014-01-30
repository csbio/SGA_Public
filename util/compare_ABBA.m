function [labels, abba_eps, abba_pvl] = compare_ABBA(sga)
%function [] = compare_ABBA(sga)

    equiv_file = '~/SGA/Main/refdata/array_query_map_equiv_good_only_140105.csv';
    fprintf('using equiv file:\n%s\n', equiv_file);
    fid = fopen(equiv_file, 'r');
    equiv = textscan(fid, '%s%s', 'Delimiter', '\t', 'ReturnOnError', false);
    fclose(fid);
    equiv{1} = lower(equiv{1}); equiv{2} = lower(equiv{2});
    equiv_map = Hash([], equiv{1}); % arrays are in col 1

    %ABBA will be an epsilon matrix, QxA, with alleles at the same index
    QQ = sga.Cannon.Orf(sga.Cannon.isQuery);
    AA = sga.Cannon.Orf(sga.Cannon.isArray);
    [QQ_head, QQ_tails] = StripOrfs(QQ);
    [AA_head, AA_tails] = StripOrfs(AA);
    CannonMap = Hash([], sga.Cannon.Orf);

    % step two match the strains by extention
    % this logic works for deletions too because they're in the file

    ts_queries = QQ_tails(ismember(QQ_tails, equiv{2}));
    ts_arrays  = AA_tails(ismember(AA_tails, equiv{1}));
    [ts_strains, ixa, ixb] = intersect(ts_queries, apply_map(equiv_map, ts_arrays, equiv{2}));
    ts_queries = ts_queries(ixa);
    ts_arrays = ts_arrays(ixb);
    % map these back to full strain ids
    ts_queries = apply_map(Hash([], QQ_tails), ts_queries, QQ);
    ts_arrays  = apply_map(Hash([], AA_tails), ts_arrays, AA);
    labels = [ts_queries ts_arrays];

    % step three, assemble the pieces:
    abba_eps = sga.eps(apply_map(CannonMap, labels(:,1)), apply_map(CannonMap, labels(:,2)));
    abba_pvl = sga.pvl(apply_map(CannonMap, labels(:,1)), apply_map(CannonMap, labels(:,2)));


% -------------------------------------------------------------------------------------------
    % Analysis workspace

    % optional: do some stats work
    % AB vs BA
    abba_eps(abba_pvl > 0.05) = NaN;
    % abba_eps(abba_eps >-0.08 & abba_eps < 0.08) = NaN;

    upper_mask = logical(triu(ones(size(abba_eps))));
    AB = abba_eps(upper_mask);
    BA = abba_eps'; % have to transpose and use the same mask b.c. col-major
    BA = BA(upper_mask);
    [r, p] = corrcoef([AB BA], 'rows', 'pairwise');
    fprintf('r = %0.2f\np = %0.2e\n', r(1,2), p(1,2));


    figure
    subplot(1,2,1);
    scatter(AB, BA, 'k.')
    SquareLine();
    xlabel('AB epsilon')
    ylabel('BA epsilon');
    subplot(1,2,2);
    hist2d(AB, BA);
    xlabel('AB epsilon')
    ylabel('BA epsilon');
    subplot(1,2,1);



end



