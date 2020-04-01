function[] = print_trigenic(sga_raw, sga_triple, filename)
%% Adapted from print_sga_2015.m and print_tau_2018.m

% alleles = StrainToAllele(sga_triple.Cannon.Orf); % May need to modify for triples

fid = fopen(filename, 'w');
header = {
    'Query Strain ID',
    %'Query allele name',
    'Array Strain ID',
    %'Array allele name',
    % 'Raw genetic interaction score (epsilon)', 
    'Adjusted genetic interaction score (epsilon or tau)',
    'P-value',
    % 'Query single mutant fitness (SMF)',
    % 'Array SMF',
    'Double/triple mutant fitness',
    'Double/triple mutant fitness standard deviation'};

h_line = header{1};
for i=2:length(header)
    h_line = [h_line sprintf('\t%s', header{i})];
end
fprintf(fid, '%s\n', h_line);

% print in row-major order
% for/find loops only work on row-vectors
% so let's check assumptions
assert(size(sga_triple.Cannon.isQuery, 2) == 1);
assert(size(sga_triple.Cannon.isArray, 1) == 1);

q_ix = find(sga_triple.Cannon.isQuery');
a_ix = find(sga_triple.Cannon.isArray);


for q=q_ix
    for a=a_ix
        % skip NaNs
        if isnan(sga_triple.eps(q,a))
            continue;
        end
        
        % first 4 text columns
        fprintf(fid, '%s\t%s\t', ...
            sga_triple.Cannon.Orf{q}, ...            
            sga_triple.Cannon.Orf{a}); 
        %alleles{q}, ...    
        %alleles{a});
        
        % score
        fprintf(fid, '%.4f\t%.3e\t', ... % sga_raw.eps(q,a), ...
            sga_triple.eps(q,a), sga_triple.pvl(q,a));
        
        % query & array fitnesses
        % fprintf(fid, '%.4f\t%.4f\t', ...
        %    sga_triple.qfit(q,a), sga_triple.afit(q,a));
        
        % double
        fprintf(fid, '%.4f\t%.4f\n', sga_triple.dbl(q,a), sga_triple.dbl_std(q,a));
    end
end
fclose(fid);

end
