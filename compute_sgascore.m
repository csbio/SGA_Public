%% 
% COMPUTE_SGASCORE - performs a series of normalizations on raw colony size data and
%   measures quantitative genetic interactions
%
% Inputs:
%	inputfile        path of raw colony size data
%	outputfile       path of output file
%	removearraylist  path of file containing the list of arrays to be removed
%	linkagefile      path of file containing the query-specific linkage window sizes
%	smfitnessfile    path of file containing single mutant fitness
%	wild_type        strain-id of wild type cols [Default = 'URA3control_sn4757']
%	border_strain_orf strain-id of border cols.   [Default = 'YOR202W_dma1']
%       eps_qnorm_ref    string to matfile containing eps_norm_table for quantile
%       random_seed	 random seed used to reproduce a random run
%
% Boolean Flags [default]:
%	skip_linkage_detection [F]: don't bother detecting linkage.
%			If this is set skip_linkage_mask has no effect either way
%
%	skip_linkage_mask [F]: if set to true, detect linkage and hold it 
%			out for corrections, but PUT IT BACK before the end.
%			Effectively, detect linkage, but don't mask it out.	
%
%	skip_perl_step [F]: The preprocessing step can be quite slow.
%			Set this to true to skip it (if you've scored the input once before)
%
%	skip_wt_remove [F]: Usually, WT-queries are NaN'ed out at the end
%			Set this to true to skip that and print out a WT query profile
%
%  remove_HO_globally [F]: The script will detect if all queries are double
%        mutants and will remove HO linkage for all queries. However, if you
%        have added other queries for stability, this detection will fail. 
%        This flag gives you the option to remove HO linkage for all queries
%        even when not all queries are DM. If this option is false, auto-detection 
%        still proceeds.
%
%  skip_batch_correction [F]: This will skip the batch correction step.
%
%  disable_arrayvar_pval [F]: This will remove the "second pvalue" derived from
%        array variance. Set to TRUE for the condition project.
%
%  disable_jackknife [F]: Disable the jackknife colony filter.
%         Set to TRUE for the condition project
%
%	
% Authors: Chad Myers (cmyers@cs.umn.edu), 
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca),
%          Benjamin VanderSluis (bvander@cs.umn.edu)
%
% Last revision: 2012-05-04
%
%% Checks before starting

compute_sgascore_tic = tic;

%% Default project head and path settings.
% only use functions from _this_ code tree
% this started throwing warnings about an optifun directory not existing
% it does exist, has not changed, and though the warning is common, i cannot reproduce it with 100% success, #maddening
% restoredefaultpath  
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
    tmp = split_by_delimiter(outputfile, '/');
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

% now that we have a log, report our base_dir
log_printf(lfid, sprintf('Working dir is %s\n', base_dir));

% Set some default flags if not defined
keyboard_confirm = false;
if(~exist('skip_linkage_detection', 'var'))
    % Skips the whole linkage procedure
    % experimental / for side projects
    log_printf(lfid, 'Using DEFAULT: skip_linkage_detection = false\n');
    skip_linkage_detection = false;
end
if(~exist('skip_linkage_mask', 'var'))
    % hold out linkage colonies for correction, then put them BACK 
    % into the outputfile. Use this to generate "nolink" data for
    % updating the linkage region definitions.
    log_printf(lfid, 'Using DEFAULT: skip_linkage_mask = false\n');
    skip_linkage_mask = false;
end
if ~exist('skip_perl_step', 'var')
    log_printf(lfid, 'Using DEFAULT: skip_perl_step = false\n');
    skip_perl_step = false;
end
if ~exist('wild_type', 'var')
    log_printf(lfid, 'Using DEFAULT: wild_type = URA3control_sn4757\n');
    wild_type = 'URA3control_sn4757';
    keyboard_confirm = true;
end
if ~exist('border_strain_orf', 'var')
    log_printf(lfid, 'Using DEFAULT: border_strain_orf = YOR202W_dma1\n');
    border_strain_orf = 'YOR202W_dma1';
    keyboard_confirm = true;
end
if ~exist('skip_wt_remove', 'var')
    log_printf(lfid, 'Using DEFAULT: skip_wt_remove = false\n');
    skip_wt_remove = false;
end
if ~exist('remove_HO_globally', 'var')
   log_printf(lfid, 'Using DEFAULT: remove_HO_globally = false\n');
   remove_HO_globally = false;
end
if ~exist('skip_batch_correction', 'var')
   log_printf(lfid, 'Using DEFAULT: skip_batch_correction = false\n');
   skip_batch_correction = false;
end
if ~exist('disable_arrayvar_pval', 'var')
   log_printf(lfid, 'Using DEFAULT: disable_arrayvar_pval = false\n');
   disable_arrayvar_pval = false;
end
if ~exist('disable_jackknife', 'var')
   log_printf(lfid, 'Using DEFAULT: disable_jackknife = false\n');
   disable_jackknife = false;
end
if ~exist('disable_big_colony_removal', 'var')
    log_printf(lfid, 'Using DEFAULT: disable_big_colony_removal = false\n');
    disable_big_colony_removal = false;
end


if(keyboard_confirm)
    log_printf(lfid, 'Automatic selections made. Please confirm.\n');
    log_printf(lfid, 'Type "dbcont" to continue or "dbquit" to abort.');
    keyboard;
end

% Check inputs
vars_to_check = {'inputfile','removearraylist','linkagefile', 'coord_file', 'smfitnessfile'};
if(skip_linkage_detection)
    vars_to_check = setdiff(vars_to_check, {'linkagefile', 'coord_file'});
end
for i = 1 : length(vars_to_check) 
    if ~exist(vars_to_check{i},'var')
        error(['Define ' vars_to_check{i} '.']);
    else
        if ~exist(eval(vars_to_check{i}), 'file')
            error(['Indicated ' vars_to_check{i} ' does not exist: ' eval(vars_to_check{i})]);
        end
    end
end

% If we are skipping the perl step, we need to check for pre-processed files as well
if(skip_perl_step)
   if ~exist([inputfile '_orfidmap'], 'file') || ~exist([inputfile '_numeric'], 'file')
      error('skip_perl_step is set to true, but _orfidmap and/or _numeric for inputfile do not exist!')
   end
end

% Set random seed if defined
if exist('random_seed','var')
    rng(random_seed);
else
   tmp = datevec(now);
   rng(tmp(6)*1000);
   random_seed = getfield(rng, 'Seed');
   clear tmp;
end

log_printf(lfid, sprintf('Random Seed is %s\n', num2str(random_seed)));

% Save the path and version for the records
pth = path;
matlab_version = version();
version_blacklist = {'2010b'}; % 2010b has a bug in svd() that screws up batch correction
if(~isempty(strfind(matlab_version, version_blacklist{1}))) % loop this if blacklist grows
    error('You are attempting to use an unsupported version of matlab %s', version_blacklist{1})
end

%% Load raw data
if(skip_perl_step)
    log_printf(lfid, 'Skipping perl preprocessing\n');
end
sgadata = load_raw_sga_data_withbatch(inputfile, skip_perl_step, lfid);
log_printf(lfid, 'Data loaded.\n');

% Border strain - SGA = YOR202W (HIS3) TSA = YMR271C (URA10)
ind_border = strmatch(border_strain_orf, sgadata.orfnames,'exact');
num_border = sum(sgadata.arrays == ind_border);
log_printf(lfid, 'using border strain %s\n', border_strain_orf);
log_printf(lfid, 'border strain array matches %d colonies (%d%%); expected (19%%)\n', ...
                 num_border, floor(100*num_border/length(sgadata.arrays)));

sgadata.short_orf_names = sgadata.orfnames;
for i=1:length(sgadata.short_orf_names)
	delim_ix = strfind(sgadata.orfnames{i}, '_');
	if(~isempty(delim_ix))
		sgadata.short_orf_names{i} = sgadata.orfnames{i}(1:delim_ix(1)-1);
	end
end

% Check for needed WT data early
wild_type_id = strmatch(wild_type, sgadata.orfnames);
if(isempty(wild_type_id))
	err_string = sprintf('\n\nTERMINAL WARNING - Cannot calculate array strain variance, no WT screens (%s) found\nWARNING\n', wild_type);
    log_printf(lfid, err_string);
    error(err_string);
end

% report the number of strains of different types from the orfmap file
% note, not all of these are mutually exclusive with all others
query_strains = sgadata.orfnames(unique(sgadata.querys));
array_strains = sgadata.orfnames(unique(sgadata.arrays));
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
log_printf(lfid, '\n\nStrain Summary:\n');
log_printf(lfid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',   'type' , strain_types{:});
log_printf(lfid, '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',   'query' , strain_type_counts(1,:));
log_printf(lfid, '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n\n', 'array' , strain_type_counts(2,:));
unique_array_plate_count = length(unique(sgadata.arrayplateids));
log_printf(lfid, 'number of unique array plates found: %d\n', unique_array_plate_count);
if ~ismember(unique_array_plate_count, [1, 4, 5, 14])
    log_printf(lfid, 'TERMINAL unknown array configuration (%d) or incorrect array plate mapping');
    return
end

%% Define colonies to ignore for some normalization steps

% Get colonies corresponding to questionable arrays
dat = importdata(removearraylist);
[t, ind1, ind2] = intersect(sgadata.orfnames, dat);
bad_array_cols = find(ismember(sgadata.arrays, ind1));
log_printf(lfid, '%d colonies ignored from "bad arrays"\n', length(bad_array_cols));

%% Speed optimization #1: Construct plateid->ind mapping

% Minimize the plateids
[all_plateids, m, n] = unique(sgadata.plateids);
all_plateids_new = (1:length(all_plateids))';
sgadata.plateids = all_plateids_new(n);

% Constructing plateid->ind map
all_plateids = unique(sgadata.plateids);
plate_id_map = cell(max(all_plateids),1);

log_printf(lfid, ['\nConstructing plateid->ind map...\n|' blanks(50) '|\n|']);
for i = 1:length(all_plateids)
    plate_id_map{all_plateids(i)} = find(sgadata.plateids == all_plateids(i));
    print_progress(lfid, length(all_plateids),i);
end
log_printf(lfid, '|\n');

%% Speed optimization #2: Construct query->ind mapping

log_printf(lfid, 'Constructing query->ind map...\n');

all_querys = unique(sgadata.querys);
query_map = cell(max(all_querys),1);

for i=1:length(all_querys),
   query_map{all_querys(i)} =  find(sgadata.querys == all_querys(i));
end

%% Speed optimization #3: Construct array->ind mapping

log_printf(lfid, 'Constructing array->ind map...\n');

all_arrays = unique(sgadata.arrays);
array_map = cell(max(all_arrays),1);

for i=1:length(all_arrays)
   array_map{all_arrays(i)}=find(sgadata.arrays == all_arrays(i));
end

%% Speed optimization #4: Create index of each group of 4 spots

sgadata.all_arrayplateids = unique(sgadata.arrayplateids);
all_arrayplateids_map = zeros(max(sgadata.all_arrayplateids),1);
all_arrayplateids_map(sgadata.all_arrayplateids) = 1:length(sgadata.all_arrayplateids);

% Map colony coordinates back to 384 format
row384 = ceil(sgadata.rows/2);
col384 = ceil(sgadata.cols/2);
ind384 = sub2ind([16 24], row384, col384);

% Generate a replicate ID unique across multiple arrayplates
sgadata.replicateid = (all_arrayplateids_map(sgadata.arrayplateids)-1)*384 + double(ind384);
sgadata.spots = sgadata.plateids*10000 + sgadata.replicateid;

%% SAFE ENTRY

% Get colonies corresponding to linkage
linkage_cols = [];
if(~skip_linkage_detection)
	linkage_cols = filter_linkage_colonies(sgadata, linkagefile, coord_file,...
		  all_querys, all_arrays, query_map, array_map, wild_type, remove_HO_globally, lfid);
	log_printf(lfid, '%d colonies identified as linkage\n', length(linkage_cols));
end
ignore_cols = unique([bad_array_cols; linkage_cols]);

%% Normalizations
% Default median colony size per plate
default_median_colsize = 510;

% Plate normalization
sgadata.colsize_platenorm = ...
    apply_plate_normalization(sgadata, 'colsize', ignore_cols, default_median_colsize, plate_id_map, lfid);

% Filter very large colonies
ind = find(sgadata.colsize_platenorm >= 1.5*default_median_colsize & ...
    sgadata.rows > 2 & sgadata.rows < 31 & ...
    sgadata.cols > 2 & sgadata.cols < 47);

all_spots = unique(sgadata.spots(ind));
num_big_colonies_per_spot = histc(sgadata.spots(ind), all_spots);
spots_to_remove = all_spots(num_big_colonies_per_spot >= 3);

sgadata.colsize_platenorm(ismember(sgadata.spots, spots_to_remove)) = NaN;

% Get colony residuals logresiduals and medians
[sgadata.residual, sgadata.logresidual, sgadata.arraymedian] = ...
    calculate_colony_residuals(sgadata, 'colsize_platenorm', plate_id_map, lfid);

% Spatial normalization
sgadata.spatialnorm_colsize = ...
    apply_spatial_normalization(sgadata, 'logresidual', ignore_cols, plate_id_map, lfid);

% Row/column correction
sgadata.rowcolcorr_colsize = ...
    apply_rowcol_normalization(sgadata, 'spatialnorm_colsize', ignore_cols, plate_id_map, lfid);

% Competition correction
sgadata.residual_spatialnorm = sgadata.rowcolcorr_colsize - sgadata.arraymedian;
sgadata.compcorr_colsize = ...
    apply_competition_correction(sgadata, 'residual_spatialnorm', ignore_cols, plate_id_map, lfid);

% One last round of plate scaling
sgadata.finalplatecorr_colsize = apply_plate_normalization(sgadata, 'compcorr_colsize', ...
                                 ignore_cols, default_median_colsize, plate_id_map, lfid);

if disable_jackknife
   sgadata.filt_jackknife = sgadata.finalplatecorr_colsize;
else
   % Do filtering based on held-out CV values
   sgadata.filt_jackknife = apply_jackknife_correction(sgadata, 'finalplatecorr_colsize', ...
                          border_strain_orf, query_map, plate_id_map, lfid);
end

% keep a copy (in jackknife) to replace linkage after batch
% SAFE ENTRY POINT, this is the first time we use filt colsize
sgadata.filt_colsize = sgadata.filt_jackknife;

%% Filters section

% Set values == 0 to NaN
sgadata.filt_colsize(sgadata.filt_colsize < 1) = NaN;

% Set values > 1000 to 1000
if (~disable_big_colony_removal)
    sgadata.filt_colsize(sgadata.filt_colsize > 1000) = 1000;
end

% Remove border (contains HIS3 control strains)
ind = find(sgadata.rows < 3 | sgadata.rows > 30 | sgadata.cols < 3 | sgadata.cols > 46);
sgadata.filt_colsize(ind) = NaN;

% Delete any array strain that is missing more than 15% of its supposed occurrences
tot_cols = zeros(length(all_arrays),1);
good_cols = zeros(length(all_arrays),1);

for i = 1:length(all_arrays)
   tot_cols(i,1) = length(array_map{all_arrays(i)});
   good_cols(i,1) = length(find(~isnan(sgadata.filt_colsize(array_map{all_arrays(i)}))));
end

ind = find(good_cols < 0.85 * tot_cols);
ind = setdiff(ind,strmatch(border_strain_orf,sgadata.orfnames(all_arrays),'exact'));

remove_ind = find(ismember(sgadata.arrays, all_arrays(ind)));
sgadata.filt_colsize(remove_ind) = NaN;


% Remove "ignore_cols" (bad arrays and linkage if applicable)
% linkage cols will be replaced later
sgadata.filt_colsize(ignore_cols) = NaN;

log_printf(lfid, 'Finished applying filters...\n');

%% Remove Triple-Specific interactions BJV
% mostly known slow-growers in uracil selection step

%URA1 URA2 URA4 URA5 URA7 URA8 URA10 EST1 MDM12 GRR1
triple_remove_list = {'YKL216W_dma2907', 'YJL130C_dma2510', 'YLR420W_dma3476', ...
                      'YML106W_dma3610', 'YBL039C_dma130',  'YJR103W_dma2710', ...
                      'YMR271C_dma3865', 'YLR233C_dma3278', 'YOL009C_dma4290', 'YJR090C_dma2675'};
	for i=1:length(triple_remove_list)
    	triple_remove_list{i} = strmatch(triple_remove_list{i}, sgadata.orfnames, 'exact');
	end
	triple_remove_list = triple_remove_list(~cellfun('isempty', triple_remove_list));
	triple_remove_list = cell2mat(triple_remove_list);
	
	% Locate the DM queries as marked with a '+'
	plus_cell = cell(size(sgadata.orfnames));
	plus_cell(:) = {'+'};
	triple_queries = find(~cellfun(@isempty, cellfun(@strfind, sgadata.orfnames, plus_cell, 'UniformOutput', false)));
	
	% Any colony which matches any combination of a query and array above is removed
	triple_queries_bool = logical(zeros(length(sgadata.querys), 1));
	triple_remove_bool  = logical(zeros(length(sgadata.arrays), 1));
	for i=1:length(triple_queries)
    	triple_queries_bool = triple_queries_bool | sgadata.querys == triple_queries(i);
	end
	for i=1:length(triple_remove_list)
    	triple_remove_bool = triple_remove_bool | sgadata.arrays == triple_remove_list(i);
	end
	
	sgadata.filt_colsize(triple_queries_bool & triple_remove_bool) = NaN;
	clear plus_cell triple_queries triple_remove_list triple_queries_bool triple_remove_bool
     

%% Batch Correction 
if skip_batch_correction
   log_printf(lfid, 'SKIPPING batch correction.\n');
   sgadata.batchnorm_colsize = sgadata.filt_colsize;
else
   sgadata = batch_correction_wrapper(sgadata, plate_id_map, lfid);
end

%% Calculate array WT variance
% We need vertcat in case wild_type_id has more than one element (i.e. replicate)
wild_type_colonies = vertcat(query_map{wild_type_id});
array_vars = zeros(length(all_arrays),2);
log_printf(lfid, ['Calculating array WT variance...\n|' blanks(50) '|\n|']);

for i = 1:length(all_arrays)
    ind = intersect(wild_type_colonies, array_map{all_arrays(i)});
    t = max(sgadata.batchnorm_colsize(ind),1);
    t(isnan(sgadata.batchnorm_colsize(ind))) = NaN;
    array_vars(i,:)=[nanmean(log(t)),nanvar(log(t))];
    print_progress(lfid, length(all_arrays),i);
end
log_printf(lfid, '|\n');

nanind = find(isnan(array_vars(:,1)));

% For the NaNs, we need to compute means from the other screens (not just WT).
% This is useful when there's a WT-screen specific issue (e.g., URA3 linkage group).

for i = nanind'
    ind = array_map{all_arrays(i)};
    t = max(sgadata.batchnorm_colsize(ind),1);
    t(isnan(sgadata.batchnorm_colsize(ind))) = NaN;
    array_vars(i,1)=nanmean(log(t));
end
array_vars(nanind,2) = nanmedian(array_vars(:,2)); % Compute vars from the median across all others

%% Calculate query variance

all_querys = unique(sgadata.querys);

% Mean and variance of double mutants
sgadata.dm_normmean = zeros(length(sgadata.colsize),1)+NaN;
sgadata.dm_normvar = zeros(length(sgadata.colsize),1)+NaN;
sgadata.dm_num = zeros(length(sgadata.colsize),1)+NaN;

sgadata.batchnorm_colsize_nonegs = max(sgadata.batchnorm_colsize, 1);
sgadata.batchnorm_colsize_nonegs(isnan(sgadata.batchnorm_colsize)) = NaN;

log_printf(lfid, ['Computing average for double mutants...\n|' blanks(50) '|\n|']);

for i = 1:length(all_querys)
    
    % BJV this matches Triple WT also
    if all_querys(i) == strmatch(strip_annotation(wild_type, 'first'),sgadata.orfnames)
        continue;
    end
    
    ind = query_map{all_querys(i)};
    ind = setdiff(ind, find(isnan(sgadata.batchnorm_colsize_nonegs)));
    
    spots = sgadata.replicateid(ind) + double(sgadata.plateids(ind))*1000;
    [all_spots, m, n] = unique(spots);
    
    % Compute average, variance and number of colonies for each group of spots
    mn = grpstats(log(sgadata.batchnorm_colsize_nonegs(ind)), spots, 'mean');
    vr = grpstats(log(sgadata.batchnorm_colsize_nonegs(ind)), spots, 'std').^2;
    nm = grpstats(sgadata.batchnorm_colsize_nonegs(ind), spots, 'numel');
    
    sgadata.dm_normmean(ind) = mn(n);
    sgadata.dm_normvar(ind) = vr(n);
    sgadata.dm_num(ind) = nm(n);
    
    % Print progress
    print_progress(lfid, length(all_querys), i);
    
end
log_printf(lfid, '|\n');

% Pool across arrayplates for each query %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAFE ENTRY
log_printf(lfid, ['Pooling across arrayplates for each query...\n|' blanks(50) '|\n|']);
query_arrplate_vars = pool_query_arrayplate_var(sgadata, all_querys, query_map, lfid);

% SAFE ENTRY

% Load the single mutant fitness file 
[file_smfit, file_smfit_std] = load_smf(smfitnessfile, sgadata, lfid);


% replace linkage colony values that got removed before batch correction
if(skip_linkage_mask)
	sgadata.batchnorm_colsize(linkage_cols) = sgadata.filt_jackknife(linkage_cols);
	log_printf(lfid, '* REPLACING linkage colonies.\n');
end


% Model fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_fits = zeros(length(all_arrays),length(sgadata.orfnames)+1) + NaN;
model_fit_std = zeros(length(all_arrays),length(sgadata.orfnames)+1) + NaN;

all_nans = find(isnan(sgadata.batchnorm_colsize));
ind_his3 = strmatch(border_strain_orf, sgadata.orfnames,'exact');

sgadata.arraymean_corrected = nan(size(sgadata.arraymedian));% saved but unused
log_printf(lfid, ['Model fitting...\n|' blanks(50) '|\n|']);
for i = 1:length(all_arrays)
    
    all_ind = array_map{all_arrays(i)};
    all_ind = setdiff(all_ind, all_nans);

    % Print progress before continue statements
    print_progress(lfid, length(all_arrays),i);
   
   
    if isempty(all_ind)
       continue;
    end
   
    if all_arrays(i) == ind_his3
       continue;
    end
     
    % sometimes we get an error here:
    % Warning: Polynomial is not unique; degree >= number of data points.
    % Generally I only get this when scoring replicate data, so I haven't much explored the ramifications. (BJV)
	 % it comes from not having any non-NaN query fitness estimates, so 0 data points to
	 % fit the trend for each array with a degree 1 polynomial
	 % resid_model becomes [0, 0]
    
    % subtract the trend between interactions and fitness
    % Query fitness is, in theory, normalized out of the plates, but there is still
    % a residual trend we want to correct for.
    %
    % For queries that have an SMF, estimate the magnitude of the effect, 
    % and reset the colony scores to the arraymean plus the residual
    this_array_colonies = sgadata.batchnorm_colsize(all_ind);
    ind2 = find(~isnan(file_smfit(sgadata.querys(all_ind))));

    col_data_to_fit = this_array_colonies(ind2); % these are pixel scores
    query_fit_for_model = file_smfit(sgadata.querys(all_ind(ind2))); % these are fitness

    % this will fit a fitness for each array. Used for 3 things:
    % 1 the adjustment below
    % 2 used with arrays that have an SMF to fit a model pixel_fit -> WT-relative
    % 3 the model from 2 will be used to scale pixel scores that don't have an SMF, and fill them in
    resid_model = polyfit(query_fit_for_model, col_data_to_fit, 1);
    col_data_prime  = resid_model(1) * query_fit_for_model + resid_model(2);
    resid = col_data_to_fit - col_data_prime;
    this_array_colonies(ind2) = nanmean(col_data_to_fit) + resid;
    
    % collapse replicate pixel scores and calculate std dev
    curr_querys = sgadata.querys(all_ind);
    uniq_querys = unique(curr_querys);
    
    % this is essentialy a colony-wise pixel epsilon
    query_effects = zeros(length(uniq_querys),1);
    query_effects_std = zeros(length(uniq_querys),1);
    
    for k = 1:length(uniq_querys)
        ind = find(curr_querys == uniq_querys(k));
        query_effects(k) = nanmean(this_array_colonies(ind));
        query_effects_std(k) = nanstd(this_array_colonies(ind));
    end
    
    % subtract the mean to convert to pixel residual
    arrmean = nanmean(query_effects);
    arrstd = nanstd(query_effects);
    query_effects = query_effects - arrmean;
   
    % put array mean in column 1
    query_effects = [arrmean; query_effects];
    query_effects_std = [arrstd; query_effects_std];
    
    % Added back (11-04-20)
    % shift to minimize the total number interactions
    res = adjust_model_fit(query_effects, query_effects_std);
    query_effects(1) = query_effects(1) + res;
    query_effects(2:end) = query_effects(2:end) - res;
   
    % put in canonical order, with array smf (pixel) in column 1
    [int, a, b] = intersect(1:length(sgadata.orfnames), uniq_querys);

    model_fits(i,a+1) = query_effects(b+1);
    model_fit_std(i,a+1) = query_effects_std(b+1);
    model_fits(i,1) = query_effects(1);
    model_fit_std(i,1) = query_effects_std(1);


    % record arraymean_corrected for posterity
    % arraymean_corrected is not used in any further calculation
    all_ind = array_map{all_arrays(i)}; % get nan_colonies back
    sgadata.arraymean_corrected(all_ind) = nanmean(query_effects + arrmean);
    
end
log_printf(lfid, '|\n');


% Fill in SM fitness for arrays that appear here that aren't in standard
model_smfits = zeros(size(file_smfit,1),1);
model_smstd = zeros(size(file_smfit,1),1);

model_smfits(all_arrays) = model_fits(:,1);
model_smstd(all_arrays) = model_fit_std(:,1);

model_smfits(model_smfits == 0) = NaN;
model_smstd(model_smstd == 0) = NaN;

% fit a model to transform pixel fitness scores scores to WT-relative
ind = find(~isnan(file_smfit(:,1)) & ~isnan(model_smfits));
p = polyfit(model_smfits(ind), file_smfit(ind),1);

ind = find(isnan(file_smfit(:,1)));
final_smfit = file_smfit;
final_smfit(ind) = model_smfits(ind)*p(1)+p(2);
final_smfit_std = file_smfit_std;
final_smfit_std(ind) = model_smstd(ind)*p(1);


% SAFE ENTRY
sga_score = zeros(length(sgadata.orfnames)) + NaN;
sga_score_std = zeros(length(sgadata.orfnames)) + NaN;
sga_score(:,all_arrays) = model_fits(:,2:end)';
sga_score_std(:,all_arrays) = model_fit_std(:,2:end)';

avar = zeros(1,length(sgadata.orfnames)) + NaN;
amean = zeros(1,length(sgadata.orfnames)) + NaN;

avar(all_arrays) = array_vars(:,2);
amean(all_arrays) = array_vars(:,1);

amat = repmat(avar,size(sga_score,1),1);
amat_mean = repmat(amean,size(sga_score,1),1);


% Create array to arrayplate map
array_arrplate = cell(max(sgadata.arrays),1);
for i = 1:length(all_arrays)
    
    ind = unique(sgadata.arrayplateids(array_map{all_arrays(i)}));
    array_arrplate{all_arrays(i)} = ind;
    
end

qmat = zeros(size(sga_score,1),size(sga_score,2)) + NaN;
for j = 1:length(all_arrays)
    
    arrplate_ind = all_arrayplateids_map(array_arrplate{all_arrays(j)}); % arrayplate(s) this array appears on
    ind = find(ismember(query_arrplate_vars(:,2), arrplate_ind));
       
    for i = 1:length(all_querys)
       ind2 = find(query_arrplate_vars(ind,1) == i);
       qmat(all_querys(i),all_arrays(j))= nanmean(query_arrplate_vars(ind(ind2),4));
    end
    
end

% upper/lower CI bound is backgd_mean*(backgd_std)^n OR backgd_mean/(backgd_std)^n where n is the number of std. dev. for normal cdf
backgd_mean = exp(amat_mean);
backgd_std = exp(sqrt(amat+qmat));  
clear amat amat_mean;



% To avoid removing real interactions, assume the query fitness is 1 if we
% don't have other fitness information. This will cause incorrectly scaled
% interactions but we have no way of dealing with this is we don't have SM
% fitness information.
qfit = repmat(final_smfit,1,length(sga_score));
qfit(isnan(qfit)) = 1;

afit = repmat(final_smfit',length(sga_score),1);

% fit another model from pixel to wt-relative space, this time using all available array data
ind = find(~isnan(model_fits(:,1)) & ~isnan(final_smfit(all_arrays)));
p = polyfit(final_smfit(all_arrays(ind)),model_fits(ind,1),1);
c = p(1);

% convert sga_score to wt-relative space, and scale by query fitness
eps = sga_score .* (qfit/c);
eps_std = sga_score_std.*(qfit/c);

% Mainly for TS data to calibrate to FG
if(exist('eps_qnorm_ref', 'var') && ~isempty(eps_qnorm_ref))
   log_printf(lfid, '!! eps_norm_ref set, beginning quantile normalization !!\n');
   eval(['load ' eps_qnorm_ref]);
   eps = quantile_normalize_from_table(eps, eps_norm_table);
end

dm_expected = afit.*qfit;
dm_actual = dm_expected + eps;
dm_actual_std = eps_std;
clear eps eps_std afit;

% calculate pvalue 1, based on replicate spot variance
% assume this variance is normal, and assess difference from zero
% based on that assumption
pval1 = normcdf(-abs(sga_score./sga_score_std));

if disable_arrayvar_pval
   pvals = pval1;
else
   % calculate pvalue 2, based on array strain variance
   % assume that array variance is log-normally distrubuted
   % then assess difference from the background mean (of that array)
   % based on that assumption
   pval2 = normcdf(-abs(log((backgd_mean + sga_score)./backgd_mean) ./ log(backgd_std)));

   % final pvalue is the geometric mean of these two
   pvals = sqrt(pval1 .* pval2);
end




% remove wild-type data
if(~skip_wt_remove)
    ind = strmatch(strip_annotation(wild_type, 'first'), sgadata.orfnames);
    sga_score(ind,:)=NaN;
    sga_score(:,ind)=NaN;
end


%% Printing out

output_interaction_data(outputfile, sgadata.orfnames, sga_score, ...
      sga_score_std, pvals, final_smfit, final_smfit_std, ...
      dm_expected, dm_actual, dm_actual_std, lfid);

%% Final save

save('-v7.3',[outputfile,'.mat']);

compute_sgascore_time = toc(compute_sgascore_tic);
log_printf(lfid, 'total time elapsed: %.2f hours\n', compute_sgascore_time/3600);
fclose(lfid);
