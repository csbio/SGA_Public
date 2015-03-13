eval(sprintf('load %s.mat', outputfile));
outputfile = [outputfile '_qvarfix'];

compute_sgascore_tic = tic;

%% Default project head and path settings.
% only use functions from _this_ code tree
restoredefaultpath  
base_dir = pwd();
addpath(base_dir)
addpath([base_dir '/IO']);
addpath([base_dir '/corrections']);
addpath([base_dir '/util']);

% Disable the pager if enabled
more off


% Check output
if ~exist('outputfile','var')
    error('Define outputfile.');
else
    tmp = split_by_delimiter('/', outputfile);
    outputdir = join_by_delimiter(tmp(1:end-1), '/');
    if(outputfile(1) == '/')
        outputdir = ['/' outputdir];
    end
    if ~exist(outputdir,'dir')
        fprintf(sprintf('\n\n!! Current dir is: %s\n', pwd()));
        fprintf(['!! Output directory ' outputdir ' does not exist.\n']);
        fprintf('!! Type "dbcont" to create it and continue or "dbquit" to abort.\n\n')
        keyboard
        mkdir(outputdir);
    end 
    try
        lfid = fopen([outputfile '.log'], 'w');
    catch
        fprintf('cannot open log file: %s\n', [outputfile '.log']);
    end
end

% SAFE ENTRY
log_printf(lfid, ['Pooling across arrayplates for each query...\n|' blanks(50) '|\n|']);
query_arrplate_vars = pool_query_arrayplate_var(sgadata, all_querys, query_map, lfid);

% SAFE ENTRY
% Load the single mutant fitness file -----------------------------------------------------------------------------------------------------
smf_fid = fopen(smfitnessfile, 'r');
fitness_data = struct();
fitness_data.raw = textscan(smf_fid, '%s%f%f', 'Delimiter', '\t', 'ReturnOnError', false);
fclose(smf_fid);
fitness_data.ORF = fitness_data.raw{1};
fitness_data.SMF = fitness_data.raw{2};
fitness_data.STD = fitness_data.raw{3};
fitness_report_header = {'Exact match', 'Partial Match', 'Not Found', 'NaN in file'};
fitness_report_counts = zeros(1,4);
fitness_hash = hash_strings(fitness_data.ORF);

final_smfit = zeros(length(sgadata.orfnames),1) + NaN;
final_smfit_std = zeros(length(sgadata.orfnames),1) + NaN;
for i=1:length(sgadata.orfnames)
    if fitness_hash.containsKey(sgadata.orfnames{i})
        final_smfit(i)     = fitness_data.SMF(fitness_hash.get(sgadata.orfnames{i}));
        final_smfit_std(i) = fitness_data.STD(fitness_hash.get(sgadata.orfnames{i}));
        fitness_report_counts(1) = fitness_report_counts(1)+1;
    elseif fitness_hash.containsKey(strip_annotation(sgadata.orfnames{i}))
        final_smfit(i)     = fitness_data.SMF(fitness_hash.get(strip_annotation(sgadata.orfnames{i})));
        final_smfit_std(i) = fitness_data.STD(fitness_hash.get(strip_annotation(sgadata.orfnames{i})));
        fitness_report_counts(2) = fitness_report_counts(2)+1;
    else
        fitness_report_counts(3) = fitness_report_counts(3)+1;
    end
end
fitness_report_counts(4) = sum(isnan(final_smfit)) - fitness_report_counts(3);

log_printf(lfid, 'Fitness file report:\n');
for i=1:length(fitness_report_header)
    log_printf(lfid, '\t%s\t: %d\n', fitness_report_header{i}, fitness_report_counts(i));
end
% Load the single mutant fitness file -----------------------------------------------------------------------------------------------------


% replace linkage colony values that got removed before batch correction
if(skip_linkage_mask)
	sgadata.batchnorm_colsize(linkage_cols) = sgadata.filt_jackknife(linkage_cols);
	log_printf(lfid, '* REPLACING linkage colonies.\n');
end

model_fits = zeros(length(all_arrays),length(sgadata.orfnames)+1) + NaN;
model_fit_std = zeros(length(all_arrays),length(sgadata.orfnames)+1) + NaN;

all_nans = find(isnan(sgadata.batchnorm_colsize));
ind_his3 = strmatch(border_strain_orf, sgadata.orfnames,'exact');

sgadata.arraymean_corrected = nan(size(sgadata.arraymedian));
log_printf(lfid, ['Model fitting...\n|' blanks(50) '|\n|']);
for i = 1:length(all_arrays)
    
    all_ind = array_map{all_arrays(i)};
    all_ind = setdiff(all_ind, all_nans);
   
    if isempty(all_ind)
       continue;
    end
   
    if all_arrays(i) == ind_his3
       continue;
    end

    querys = unique(sgadata.querys(all_ind));
     
    ind2 = find(~isnan(final_smfit(sgadata.querys(all_ind))));
    % sometimes we get an error here:
    % Warning: Polynomial is not unique; degree >= number of data points.
    % Generally I only get this when scoring replicate data, so I haven't much explored the ramifications. (BJV)
	% it comes from not having any non-NaN query fitness estimates, so 0 data points to
	% fit the trend for each array with a degree 1 polynomial
	% p becomes [0, 0]
    p = polyfit(final_smfit(sgadata.querys(all_ind(ind2))),sgadata.batchnorm_colsize(all_ind(ind2)),1);
    
    curr_data = sgadata.batchnorm_colsize(all_ind);
    
    % Added back (11-03-31)
    % subtract the trend between interactions and fitness
    curr_data(ind2) = curr_data(ind2) + (nanmean(sgadata.batchnorm_colsize(all_ind(ind2))) - (p(1)*final_smfit(sgadata.querys(all_ind(ind2)))+p(2)));

    curr_querys = sgadata.querys(all_ind);
    uniq_querys = unique(curr_querys);
    
    query_effects = zeros(length(uniq_querys),1);
    query_effects_std = zeros(length(uniq_querys),1);
    
    for k = 1:length(uniq_querys)
        
        ind = find(curr_querys == uniq_querys(k));
        query_effects(k) = nanmean(curr_data(ind));
        query_effects_std(k) = nanstd(curr_data(ind));
    end
    
    
    arrmean = nanmean(query_effects);
    arrstd = nanstd(query_effects);
    query_effects = query_effects - arrmean;
   
    query_effects = [arrmean; query_effects];
    query_effects_std = [arrstd; query_effects_std];
    
    % Added back (11-04-20)
    % shift to minimize the total number interactions
    res = adjust_model_fit(query_effects, query_effects_std);
    query_effects(1) = query_effects(1) + res;
    query_effects(2:end) = query_effects(2:end) - res;
   
    [int, a, b] = intersect(1:length(sgadata.orfnames), uniq_querys);

    model_fits(i,a+1) = query_effects(b+1);
    model_fit_std(i,a+1) = query_effects_std(b+1);
    model_fits(i,1) = query_effects(1);
    model_fit_std(i,1) = query_effects_std(1);

    all_ind = array_map{all_arrays(i)}; % get nan_colonies back to record corrected mean
    % put the removed mean back in and add the query_effects???
    sgadata.arraymean_corrected(all_ind) = nanmean(query_effects + arrmean);
    
    % Print progress
    print_progress(lfid, length(all_arrays),i);
   
end
log_printf(lfid, '|\n');


% Fill in SM fitness for arrays that appear here that aren't in standard
model_smfits = zeros(size(final_smfit,1),1);
model_smstd = zeros(size(final_smfit,1),1);

model_smfits(all_arrays)=model_fits(:,1);
model_smstd(all_arrays)=model_fit_std(:,1);

model_smfits(model_smfits == 0) = NaN;
model_smstd(model_smstd == 0) = NaN;

ind = find(~isnan(final_smfit(:,1)) & ~isnan(model_smfits));
p = polyfit(model_smfits(ind), final_smfit(ind),1);

ind = find(isnan(final_smfit(:,1)));
final_smfit(ind) = model_smfits(ind)*p(1)+p(2);
final_smfit_std(ind) = model_smstd(ind)*p(1);


% SGA score
% SAFE ENTRY
sga_score = model_fits(:,2:end);
sga_score_std = model_fit_std(:,2:end);

complete_mat = zeros(length(sgadata.orfnames)) + NaN;
complete_mat_std = zeros(length(sgadata.orfnames)) + NaN;
complete_mat(:,all_arrays) = sga_score';
complete_mat_std(:,all_arrays) = sga_score_std';
clear sga_score sga_score_std;

avar = zeros(1,length(sgadata.orfnames)) + NaN;
amean = zeros(1,length(sgadata.orfnames)) + NaN;

avar(all_arrays) = array_vars(:,2);
amean(all_arrays) = array_vars(:,1);

amat = repmat(avar,size(complete_mat,1),1);
amat_mean = repmat(amean,size(complete_mat,1),1);


% Create array to arrayplate map
array_arrplate = cell(max(sgadata.arrays),1);
for i = 1:length(all_arrays)
    
    ind = unique(sgadata.arrayplateids(array_map{all_arrays(i)}));
    array_arrplate{all_arrays(i)} = ind;
    
end

qmat = zeros(size(complete_mat,1),size(complete_mat,2)) + NaN;
for j = 1:length(all_arrays)
    
    arrplate_ind = all_arrayplateids_map(array_arrplate{all_arrays(j)}); % arrayplate(s) this array appears on
    ind = find(ismember(query_arrplate_vars(:,2), arrplate_ind));
       
    for i = 1:length(all_querys)
       ind2 = find(query_arrplate_vars(ind,1) == i);
       qmat(all_querys(i),all_arrays(j))= nanmean(query_arrplate_vars(ind(ind2),4));
    end
    
end

backgd_mean = exp(amat_mean);
backgd_std = exp(sqrt(amat+qmat));  % upper/lower CI bound is backgd_mean*(backgd_std)^n OR backgd_mean/(backgd_std)^n where n is the number of std. dev. for normal cdf
clear amat amat_mean;


qfit = repmat(final_smfit,1,length(complete_mat));

% To avoid removing real interactions, assume the query fitness is 1 if we
% don't have other fitness information. This will cause incorrectly scaled
% interactions but we have no way of dealing with this is we don't have SM
% fitness information.

qfit(isnan(qfit)) = 1;

afit = repmat(final_smfit',length(complete_mat),1);

dm_expected = afit.*qfit;
ind = find(~isnan(model_fits(:,1)) & ~isnan(final_smfit(all_arrays)));
p = polyfit(final_smfit(all_arrays(ind)),model_fits(ind,1),1);
c = p(1);

eps = complete_mat .* (qfit/c);
eps_std = complete_mat_std.*(qfit/c);

% Mainly for TS data to calibrate to FG
if(exist('eps_qnorm_ref', 'var') && ~isempty(eps_qnorm_ref))
   log_printf(lfid, '!! eps_norm_ref set, beginning quantile normalization !!\n');
   eval(['load ' eps_qnorm_ref]);
   eps = quantile_normalize_from_table(eps, eps_norm_table);
end

dm_actual = dm_expected + eps;
dm_actual_std = eps_std;
clear eps eps_std afit;

% remove wild-type data
if(~skip_wt_remove)
    ind = strmatch(strip_annotation(wild_type, 'first'), sgadata.orfnames);
    complete_mat(ind,:)=NaN;
    complete_mat(:,ind)=NaN;
end


%% Printing out

output_interaction_data(outputfile,sgadata.orfnames,complete_mat,complete_mat_std,backgd_mean,backgd_std,final_smfit,final_smfit_std,dm_expected,dm_actual,dm_actual_std,lfid);

%% Final save

save('-v7.3',[outputfile,'.mat']);

compute_sgascore_time = toc(compute_sgascore_tic);
log_printf(lfid, 'total time elapsed: %.2f hours\n', compute_sgascore_time/3600);
fclose(lfid);
