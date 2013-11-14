function[stats] = BasicSgaStats(sga, thresh, holdout)
%function[stats] = BasicSgaStats(sga, thresh, holdout)
% stats = pos_int neg_int +/-

if(~exist('holdout', 'var'))
	holdout = [];
end
if(~exist('thresh', 'var'))
	thresh = 0.08;
end

sga.eps(holdout,:) = [];
sga.eps(:,holdout) = [];
sga.pvl(holdout,:) = [];
sga.pvl(:,holdout) = [];
sga.Cannon.isQuery(holdout) = [];
sga.Cannon.isArray(holdout) = [];


EPS = sga.eps(sga.Cannon.isQuery, sga.Cannon.isArray);
PVL = sga.pvl(sga.Cannon.isQuery, sga.Cannon.isArray);
SCR = ~isnan(EPS) & ~isnan(PVL); 

num_nans = sum(sum(isnan(EPS)));

EPS(PVL > 0.05) = 0;
EPS(EPS > -thresh & EPS < thresh) = 0;

num_neg = sum(sum(EPS < 0));
num_pos = sum(sum(EPS > 0));
num_total = num_neg+num_pos;
neg_ratio = num_neg / num_total;
stats = [num_pos, num_neg, num_pos / num_neg];

%neg_den = num_neg / (sum(sga.Cannon.isQuery) * sum(sga.Cannon.isArray));
%pos_den = num_pos / (sum(sga.Cannon.isQuery) * sum(sga.Cannon.isArray));
neg_den = num_neg / sum(sum(SCR));
pos_den = num_pos / sum(sum(SCR));

fprintf('-----------------\n');
fprintf('%d queries\n%d arrays\n', sum(sga.Cannon.isQuery), sum(sga.Cannon.isArray));
fprintf('%d interactions screened\n', sum(sum(SCR)));
fprintf('interactions (- %d)(+ %d)\n', num_neg, num_pos);
fprintf('-den: %.2f%%    +den: %.2f%%\n', neg_den*100, pos_den*100)
fprintf('%.2f ratio (neg / total)\n', neg_ratio);



fprintf('--\n\n');

query_strains = sga.Cannon.Orf(sga.Cannon.isQuery);
array_strains = sga.Cannon.Orf(sga.Cannon.isArray);
strain_types = {'sn' 'dma' 'tsq' 'damp' 'tsa' 'trip' 'unann' 'total'};
strain_type_counts = nan(2,8); % query,array ; type
    strain_type_counts(1,1) = sum(~cellfun(@isempty, strfind(query_strains, '_sn')));
    strain_type_counts(1,2) = sum(~cellfun(@isempty, strfind(query_strains, '_dma')));
    strain_type_counts(1,3) = sum(~cellfun(@isempty, strfind(query_strains, '_tsq')));
    strain_type_counts(1,4) = sum(~cellfun(@isempty, strfind(query_strains, '_damp')));
    strain_type_counts(1,5) = sum(~cellfun(@isempty, strfind(query_strains, '_tsa')));
    strain_type_counts(1,6) = sum(~cellfun(@isempty, strfind(query_strains, '+')));
    strain_type_counts(1,7) = sum(cellfun(@isempty, strfind(query_strains, '_')));
    strain_type_counts(1,8) = length(query_strains);

    strain_type_counts(2,1) = sum(~cellfun(@isempty, strfind(array_strains, '_sn')));
    strain_type_counts(2,2) = sum(~cellfun(@isempty, strfind(array_strains, '_dma')));
    strain_type_counts(2,3) = sum(~cellfun(@isempty, strfind(array_strains, '_tsq')));
    strain_type_counts(2,4) = sum(~cellfun(@isempty, strfind(array_strains, '_damp')));
    strain_type_counts(2,5) = sum(~cellfun(@isempty, strfind(array_strains, '_tsa')));
    strain_type_counts(2,6) = sum(~cellfun(@isempty, strfind(array_strains, '+')));
    strain_type_counts(2,7) = sum(cellfun(@isempty, strfind(array_strains, '_')));
    strain_type_counts(2,8) = length(array_strains);
fprintf('Strain Summary:\n');
fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',   'type' , strain_types{:});
fprintf('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',   'query' , strain_type_counts(1,:));
fprintf('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n\n', 'array' , strain_type_counts(2,:));
fprintf('-----------------------------------------------------------\n');
