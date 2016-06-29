function[stats] = BasicSgaStats(sga, thresh)
%function[stats] = BasicSgaStats(sga, thresh)
% stats = pos_int neg_int +/-
% stats struct:
% .pos  positive interaction count
% .neg   negative itneraction count
% .scr  screened interactions (non-nan)
% .ratio  neg / total

if(~exist('thresh', 'var'))
	thresh = 0.08;
end

EPS = sga.eps(sga.Cannon.isQuery, sga.Cannon.isArray);
PVL = sga.pvl(sga.Cannon.isQuery, sga.Cannon.isArray);
SCR = ~isnan(EPS) & ~isnan(PVL); 

num_neg = sum(sum(EPS <= -thresh & PVL < 0.05));
num_pos = sum(sum(EPS >= thresh & PVL < 0.05));
num_total = num_neg+num_pos;
neg_ratio = num_neg / num_total;

stats = struct();
stats.pos = num_pos;
stats.neg = num_neg;
stats.scr = sum(sum(SCR));
stats.ratio = neg_ratio;
stats.queries = struct();
stats.arrays  = struct();

stats.neg_den = num_neg / sum(sum(SCR));
stats.pos_den = num_pos / sum(sum(SCR));

fprintf('-----------------\n');
fprintf('%d queries\n%d arrays\n', sum(sga.Cannon.isQuery), sum(sga.Cannon.isArray));
fprintf('%d interactions screened\n', sum(sum(SCR)));
fprintf('interactions (- %d)(+ %d)\n', num_neg, num_pos);
fprintf('-den: %.2f%%    +den: %.2f%%\n', stats.neg_den*100, stats.pos_den*100)
fprintf('%.2f ratio (neg / total)\n', neg_ratio);
fprintf('----\n\n');

query_strains = sga.Cannon.Orf(sga.Cannon.isQuery);
query_strains_common = sga.Cannon.Common(sga.Cannon.isQuery);
array_strains = sga.Cannon.Orf(sga.Cannon.isArray);
strain_types = {'sn' 'dma' 'tsq' 'damp' 'tsa' 'trip' 'y' 'supp' 'no_' 'total'};
strain_ident = {'_sn' '_dma' '_tsq' '_damp' '_tsa' '+' '_y' '_S'};

strain_type_counts = nan(2,length(strain_types)); % query,array ; type

for i=1:length(strain_ident)
    strain_type_counts(1,i) = sum(~cellfun(@isempty, strfind(query_strains, strain_ident{i})));
    strain_type_counts(2,i) = sum(~cellfun(@isempty, strfind(array_strains, strain_ident{i})));
end

% these two are custom
strain_type_counts(1,end-1) = sum(cellfun(@isempty, strfind(query_strains, '_')));
strain_type_counts(1,end) = length(query_strains);

strain_type_counts(2,end-1) = sum(cellfun(@isempty, strfind(array_strains, '_')));
strain_type_counts(2,end) = length(array_strains);

% put the results in the struct
for i=1:length(strain_types)
    stats.queries.(strain_types{i}) = strain_type_counts(1,i);
    stats.arrays.(strain_types{i})  = strain_type_counts(2,i);
end

% print them out pretty
fprintf('Strain Summary:\n');
print_table = [[{'type'}, strain_types]; [{'query', 'array'}' num2cell(strain_type_counts)]];
disp(print_table);

fprintf('-----------------------------------------------------------\n');

% list collab query strains by common name
fprintf(1, 'Collab Queries (_y):\n');
collabs = query_strains_common(substrmatch('_y', query_strains_common));
collabs = [collabs; cell(mod(4-mod(length(collabs),4),4),1)];
collabs = reshape(collabs, length(collabs)/4, 4);
cell2csv(1, collabs);



