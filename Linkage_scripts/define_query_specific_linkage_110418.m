% Input: 
%   dataset.queries     Cell array of query ORFs (with or without annot.)
%   dataset.arrays      Cell array of array ORFs (with or without annot.)
%   dataset.scores_eps  Matrix containing the genetic interaction scores
%
% Output:
%   lnkg.orf            Cell array of query ORFs (as in input)
%   lnkg.coord_mean          Matrix (queries x 2) with the linkage boundaries

function lnkg = define_query_specific_linkage_110418(dataset)


%% Load the data
labels_row = dataset.queries;
labels_col = dataset.arrays;
data = dataset.scores_eps;
[labels_row, ix1] = sort_by_chr(labels_row);
[labels_col, ix2] = sort_by_chr(labels_col);
data = data(ix1, ix2);

chromosomes = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'};

screens = labels_row;
[labels_row, annotation] = strtok(labels_row,'_');

arrays = labels_col;
[labels_col, annotation] = strtok(labels_col,'_');

load chr_length_110207;

load orf_coordinates_110207;
coord = orf_coord; clear orf_coord;

%% Get the chromosomal coordinates of queries and arrays
inds1 = [1 : length(labels_row)]';
inds2 = multistrmatch(labels_row, coord.orf,1,1,1);
ii = find(inds2 == 0);
inds1(ii) = []; 
inds2(ii) = [];

labels_row = labels_row(inds1);
data = data(inds1,:);
coord_row = [coord.start(inds2) coord.end(inds2)];

inds1 = [1 : length(labels_col)]';
inds2 = multistrmatch(labels_col, coord.orf,1,1,1);
ii = find(inds2 == 0);
inds1(ii) = [];
inds2(ii) = [];

labels_col = labels_col(inds1);
data = data(:,inds1);
coord_col = [coord.start(inds2) coord.end(inds2)];

lnkg.orf = screens;
lnkg.coord_mean = zeros(length(screens),2)+NaN;
%lnkg.coord_median = zeros(length(labels_row),2)+NaN;

%%
for ic = 1:16 
    
    iq = strmatch(['Y' chromosomes(ic)], labels_row);
    ia = strmatch(['Y' chromosomes(ic)], labels_col);
    
    all_data(ic).labels_row = labels_row(iq);
    all_data(ic).labels_col = labels_col(ia);
    all_data(ic).coord_row = coord_row(iq,:);
    all_data(ic).coord_col = coord_col(ia,:);
    all_data(ic).scores = data(iq,ia);
    all_data(ic).coord_col_length = coord_col(ia,2)-coord_col(ia,1)+1;
    
    linkage_profile = zeros(length(iq), chr_length.length(ic));
    
    % The size of the window being averaged (in bp)
    smooth_window2 = 40000;
    
    % Quite complicated way to speed up the smooth process
    coverage.bp = zeros(length(ia),chr_length.length(ic));
    
    for a = 1 : length(ia)
        
        st = all_data(ic).coord_col(a,1) - smooth_window2/2 + 1;
        en = all_data(ic).coord_col(a,2) + smooth_window2/2;
        
        coord_partial_left = st+1:st+all_data(ic).coord_col_length(a)-1;
        coord_partial_right = en-(all_data(ic).coord_col_length(a)-1):en-1;
        coord_complete = max(coord_partial_left)+1:min(coord_partial_right)-1;
       
        partial_left = 1:all_data(ic).coord_col_length(a)-1;
        partial_right = all_data(ic).coord_col_length(a)-1:-1:1;
        complete = zeros(1,length(coord_complete))+all_data(ic).coord_col_length(a);
                
        coverage.bp(a,coord_partial_left(coord_partial_left>0)) = partial_left(coord_partial_left>0);
        coverage.bp(a,coord_partial_right(coord_partial_right<=chr_length.length(ic))) = partial_right(coord_partial_right<=chr_length.length(ic));
        coverage.bp(a,coord_complete(coord_complete>0 & coord_complete<=chr_length.length(ic))) = complete(coord_complete>0 & coord_complete<=chr_length.length(ic));
       
    end
    
    coverage.bp_binary = coverage.bp > 0;
    coverage.bp_sum = sum(coverage.bp,1);
%     coverage.bp_midpoint = ceil(coverage.bp_sum/2);
%     coverage.bp_midpoint_mat = repmat(coverage.bp_midpoint, length(ia),1);
    
%     all_data(ic).scores_binary = zeros(size(all_data(ic).scores)) + NaN;
%     inds_pos = find(all_data(ic).scores>0);
%     all_data(ic).scores_binary(inds_pos) = 1;
%     inds_neg = find(all_data(ic).scores<0);
%     all_data(ic).scores_binary(inds_neg) = -1;
%     inds_zero = find(all_data(ic).scores==0);
%     all_data(ic).scores_binary(inds_zero) = 0;

fprintf(['\nChr ' num2str(ic) ' - Function in progress...\n|' blanks(50) '|\n|']);
    
    for q = 1 : length(iq)
        % Calculating the mean
        ix = find(all_data(ic).scores(q,:) > 0.08); all_data(ic).scores(q,ix) = 0.08;   % eliminate the most positive scores to reduce the noise
        scores = repmat(all_data(ic).scores(q,:)',1,chr_length.length(ic));        
        linkage_profile(q,:) = sum(coverage.bp .* scores,1)./coverage.bp_sum;
        
%         % Calculating the median (2) - It works but the mean seems to be
%         % better
%         scores = repmat(all_data(ic).scores_binary(q,:)',1,chr_length.length(ic));
%         profile(q).median = 100 .* sum(coverage.bp .* scores,1)./coverage.bp_sum;
        
%         % Calculating the median (1) - Too long, (2) gets the same answer
%         faster
%         profile(q).median = zeros(1, chr_length.length(ic))+NaN;
%         [t,ix] = sort(all_data(ic).scores(q,:)','ascend');
%         
%         coverage.bp_ranked = coverage.bp(ix,:);
%         coverage.bp_ranked_cumsum = cumsum(coverage.bp_ranked);
%         
%         tmp = coverage.bp_ranked_cumsum - coverage.bp_midpoint_mat;
%         tmp(tmp<0) = NaN;
%         [mn,mnix] = min(tmp,[],1);
%         profile(q).median = all_data(ic).scores(q,ix(mnix));
%         profile(q).median(coverage.bp_sum == 0) = NaN;  % mask positions with zero coverage
        
        
        % Calculating the linkage window boundaries
        %boundary_left = max(find(profile(q).median(1:all_data(ic).coord_row(q,1)) > 0));
        boundary_left2 = max(find(linkage_profile(q,1:all_data(ic).coord_row(q,1)) > 0));

        if isempty(boundary_left2)
            boundary_left2 = 1;
        end
        
%        boundary_right = all_data(ic).coord_row(q,2) + min(find(profile(q).median(all_data(ic).coord_row(q,2):end)>0));
        boundary_right2 = all_data(ic).coord_row(q,2) + min(find(linkage_profile(q, all_data(ic).coord_row(q,2):end)>0));

        if isempty(boundary_right2)
            boundary_right2 = chr_length.length(ic);
        end
        
%% Figure plotting. If uncommented, will create a plot for each query. Good to try out the first time you run a linkage analysis on a new dataset.
%         figure(q);
%         hold on;
%         for a = 1 : length(ia)
%             patch([all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,2) all_data(ic).coord_col(a,2)],[0 all_data(ic).scores(q,a) all_data(ic).scores(q,a) 0],'k');
%         end
%         plot(1:chr_length.length(ic), linkage_profile(q,:),'r-','LineWidth',3);
% %        plot(1:chr_length.length(ic), profile(q).median,'c-','LineWidth',3);
%         
%         plot([boundary_left2 boundary_left2], [-100 50],'r-','LineWidth',3); 
%         plot([boundary_right2 boundary_right2], [-100 50],'r-','LineWidth',3);
%         
% 
%         %plot([chr_length.length(ic)-smooth_window2 chr_length.length(ic)],[-200 -200], 'k-','LineWidth',2);
%         %plot([boundary_right-smooth_window2/2 boundary_right+smooth_window2/2],[-100 -100],'k-','LineWidth',2);
%         grid on;
%         set(gca,'XLim',[1 chr_length.length(ic)]);
%         
%         title(all_data(ic).labels_row(q));
%         linkage_profile(q, 1:boundary_left2-1) = 0;
%         linkage_profile(q, boundary_right2+1:end) = 0;
% %         
          
        lnkg.coord_mean(iq(q),:) = [boundary_left2 boundary_right2];
%        lnkg.coord_median(iq(q),:) = [boundary_left boundary_right];
        
        print_progress(length(iq), q);
        
    end
    
    
%     query_linkage_profile.query = labels_row(iq);
%     query_linkage_profile.profile_sparse = sparse(linkage_profile);
%     clear linkage_profile;
% 
%     save(['Laboratory/SGA/Recombination_map/11-09-23/query_linkage_profile_' num2str(ic) '_sparse.mat'],'query_linkage_profile','-v7.3');
%     clear query_linkage_profile;
    
    fprintf('|\n');
end

    