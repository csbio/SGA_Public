%% 
% COMPUTE_SGASCORE - performs a series of normalizations on raw colony size data and
%   measures quantitative genetic interactions
%
% Requirements: Statistics Toolbox
%
% Inputs:
%   inputfile - path to the raw colony size data
%   outputfile - path to the output file
%   removearraylist - path to the file containing the list of arrays to be removed
%   linkagefile - path to the file containing the query-specific linkage window sizes
%   smfitnessfile - path to the file containing single mutant fitness
%
% Assumptions:
% 1) WT control screens queries: "undefined" and/or "undefined+YDL227C"
% 2) any query can have an annotation after the underscore sign ("YAL002C_test")
% 3) queries and arrays are yeast ORFs (they will get mapped to specific chromosomes and chromosomal coordinates)
% 4) the border of the plates are composed of a control WT strain (HIS3 = YOR202W) BJV See line ~80
% 5) all plates are in 1536 format (32 rows x 48 columns)
% 6) queries containing a + are assumed to be double mutants
%
% Authors: Chad Myers (cmyers@cs.umn.edu), 
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca),
%          Benjamin VanderSluis (bvander@cs.umn.edu)
%
% Last revision: 2011-12-06
%
%% Checks before starting

compute_sgascore_tic = tic;
% Check inputs
vars_to_check = {'inputfile','removearraylist','linkagefile','smfitnessfile'};
for i = 1 : length(vars_to_check) 
    if ~exist(vars_to_check{i},'var')
        error(['Define' vars_to_check{i} '.']);
    else
        if ~exist(eval(vars_to_check{i}), 'file')
            error(['Indicated ' vars_to_check{i} ' does not exist: ' eval(vars_to_check{i})]);
        end
    end
end

if(~exist('skip_linkage_mask', 'var'))
    skip_linkage_mask = false;
end
if ~exist('skip_perl_step', 'var')
    skip_perl_step = false;
end

% Check output
if ~exist('outputfile','var')
    error('Define outputfile.');
else
    tmp = split_by_delimiter('/', outputfile);
    outputdir = join_by_delimiter(tmp(1:end-1), '/');
    if ~exist(outputdir,'dir')
        error(['Output directory ' outputdir ' does not exist.']);
    else
        try
            lfid = fopen([outputfile '.log'], 'w');
        catch
            fprintf('cannot open log file: %s\n', [outputfile '.log']);
        end
    end
end

% Save the path and version for the records
pth = path;
matlab_version = version();
version_blacklist = {'2010b'};
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
if ~exist('border_strain_orf', 'var')
    border_strain_orf = 'YOR202W_dma1';
end

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
if ~ismember(unique_array_plate_count, [4, 14])
    log_printf(lfid, 'TERMINAL unknown array configuration (%d) or incorrect array plate mapping');
    return
end

% Seed the random number generator % may throw error in older matlab versions
RSTREAM = RandStream.create('mt19937ar','Seed', 'shuffle');
RandStream.setGlobalStream(RSTREAM); clear RSTREM;

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

all_arrayplateids = unique(sgadata.arrayplateids);
all_arrayplateids_map = zeros(max(all_arrayplateids),1);
all_arrayplateids_map(all_arrayplateids) = 1:length(all_arrayplateids);

% Map colony coordinates back to 384 format
row384 = ceil(sgadata.rows/2);
col384 = ceil(sgadata.cols/2);
ind384 = sub2ind([16 24], row384, col384);

% Generate a replicate ID unique across multiple arrayplates
sgadata.replicateid = (all_arrayplateids_map(sgadata.arrayplateids)-1)*384 + double(ind384);
sgadata.spots = sgadata.plateids*10000 + sgadata.replicateid;

% Get colonies corresponding to linkage
linkage_cols = filter_all_linkage_colonies_queryspecific(sgadata, linkagefile, ...
     all_querys, all_arrays, query_map, array_map, lfid);
log_printf(lfid, '%d colonies identified as linkage\n', length(linkage_cols));

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

% Do filtering based on held-out CV values
sgadata.filt_jackknife = apply_jackknife_correction(sgadata, 'finalplatecorr_colsize', ...
                       border_strain_orf, query_map, plate_id_map, lfid);

% keep a copy (in jackknife) to replace linkage after batch
sgadata.filt_colsize = sgadata.filt_jackknife;


%% Filters section

field = 'filt_colsize';

% Set values == 0 to NaN
sgadata.(field)(sgadata.(field) < 1) = NaN;

% Set values > 1000 to 1000
sgadata.(field)(sgadata.(field) > 1000) = 1000;

% Remove border (contains HIS3 control strains)
ind = find(sgadata.rows < 3 | sgadata.rows > 30 | sgadata.cols < 3 | sgadata.cols > 46);
sgadata.(field)(ind) = NaN;

% Delete any array strain that is missing more than 15% of its supposed occurrences
all_arrays = unique(sgadata.arrays);

tot_cols = zeros(length(all_arrays),1);
good_cols = zeros(length(all_arrays),1);

for i = 1:length(all_arrays)
   tot_cols(i,1) = length(array_map{all_arrays(i)});
   good_cols(i,1) = length(find(~isnan(sgadata.(field)(array_map{all_arrays(i)}))));
end

ind = find(good_cols < 0.85 * tot_cols);
ind = setdiff(ind,strmatch(border_strain_orf,sgadata.orfnames(all_arrays),'exact'));

remove_ind = find(ismember(sgadata.arrays, all_arrays(ind)));
sgadata.(field)(remove_ind) = NaN;


% Remove "ignore_cols" (bad arrays and linkage if applicable)
% linkage cols will be replaced later
sgadata.(field)(ignore_cols) = NaN;

log_printf(lfid, 'Finished applying filters...\n');

%% Remove Triple-Specific interactions BJV
% mostly known slow-growers in uracil selection step

%URA1 URA2 URA4 URA5 URA7 URA8 URA10 EST1 MDM12 GRR1
triple_remove_list = {'YKL216W', 'YJL130C',  'YLR420W',  'YML106W',  'YBL039C',  ...
                     'YJR103W',  'YMR271C',  'YLR233C',  'YOL009C',  'YJR090C'};
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
	triple_queries_bool = boolean(zeros(length(sgadata.querys), 1));
	triple_remove_bool  = boolean(zeros(length(sgadata.arrays), 1));
	for i=1:length(triple_queries)
    	triple_queries_bool = triple_queries_bool | sgadata.querys == triple_queries(i);
	end
	for i=1:length(triple_remove_list)
    	triple_remove_bool = triple_remove_bool | sgadata.arrays == triple_remove_list(i);
	end
	
	sgadata.(field)(triple_queries_bool & triple_remove_bool) = NaN;
	clear plus_cell triple_queries triple_remove_list triple_queries_bool triple_remove_bool
     

%% Batch correction
% Get array plate means (use all screens, including WT screens)

all_arrplates = unique(sgadata.arrayplateids);
arrplate_ind = zeros(max(all_arrplates),1);
arrplate_ind(all_arrplates) = 1:length(all_arrplates);  

width = 48;
height = 32;

array_means = zeros(length(all_arrplates),width*height)+NaN;

log_printf(lfid, ['Getting arrayplate means...\n|' blanks(50) '|\n|']);
for i = 1:length(all_arrplates)
    
    currplates = unique(sgadata.plateids(sgadata.arrayplateids == all_arrplates(i)));
    t = [];
    
    for j = 1:length(currplates)
        
        ind = plate_id_map{currplates(j)};
        iii = sub2ind([32,48], sgadata.rows(ind), sgadata.cols(ind));
        
        d = zeros(32,48);
        d(:,[1 2 47 48]) = NaN;
        d([1 2 31 32],:) = NaN;
        d(iii) = sgadata.(field)(ind);
        
        t = [t,d(:)];
        
    end
    
    array_means(i,:) = nanmean(t,2);
    
    % Print progress
    print_progress(lfid, length(all_arrplates),i);
    
end
log_printf(lfid, '|\n');
    

% Do batch normalization using LDA

save_mats=struct;

log_printf(lfid, ['Preparing for batch normalization...\n|' blanks(50) '|\n|']);
for i=1:length(all_arrplates)
    
    curr_plates = unique(sgadata.plateids(sgadata.arrayplateids == all_arrplates(i)));
    
    curr_mat = zeros(length(curr_plates),width*height)+NaN;
    curr_ind_mat = zeros(length(curr_plates),width*height)+NaN;
    curr_batch = zeros(length(curr_plates),1);
    
    for j = 1:length(curr_plates) 
        
        ind = plate_id_map{curr_plates(j)};
        iii = sub2ind([32,48], sgadata.rows(ind), sgadata.cols(ind));
        
        d = zeros(32,48);
        d(:,[1 2 47 48]) = NaN;
        d([1 2 31 32],:) = NaN;
        d_map = zeros(32,48) + 1; % default reference to the 1st index
        
        d_map(iii) = ind;
        d(iii) = sgadata.(field)(ind);
        
        curr_mat(j,:)=d(:);
        curr_ind_mat(j,:)=d_map(:);
        
        curr_batch(j) = unique(sgadata.batch(ind));
        
    end
    
    t = curr_mat - repmat(array_means(i,:),size(curr_mat,1),1);
           
    save_mats(i).mat = t;
    save_mats(i).mat_ind = curr_ind_mat;
    save_mats(i).batch = curr_batch;
    
    % Print progress
    print_progress(lfid, length(all_arrplates),i);
    
end
log_printf(lfid, '|\n');

    
sgadata.batchnorm_colsize = sgadata.(field);
outfield = 'batchnorm_colsize';

% Normalize out batch effect. Method: LDA (supervised)
% MRECK BATCH COLLAPSE, remove when done, burn after reading
MRECK_BATCH = max(sgadata.batch)+1;
MRECK_QUERIES = find(~cellfun(@isempty, strfind(sgadata.orfnames, '_collab')));
for i=1:length(MRECK_QUERIES)
	sgadata.batch(query_map{all_querys(MRECK_QUERIES(i))}) = MRECK_BATCH;
end
log_printf(lfid, 'Re-batching MRECK queryeies (anything marked collab):\n');
log_printf(lfid, '%d _collab queries\n\n', length(MRECK_QUERIES));

perc_var = 0.1;

log_printf(lfid, ['Batch normalization...\n|', blanks(50), '|\n|']);

for j = 1:length(all_arrplates)
    
    t = save_mats(j).mat;
    t(isnan(t)) = 0;
    curr_batch = save_mats(j).batch; 
  
    % Check if some of the batches are too small -- make a larger orphan batch
    all_batches = unique(curr_batch);
    num_batch = histc(curr_batch, all_batches);
    orphan_batch = max(curr_batch)+1;
    curr_batch(ismember(curr_batch,all_batches(num_batch < 3))) = orphan_batch;
  
    tnorm = multi_class_lda(t,curr_batch,perc_var);
    sgadata.(outfield)(save_mats(j).mat_ind(:)) = sgadata.(outfield)(save_mats(j).mat_ind(:)) + (tnorm(:)-t(:));
    
    % Print progress
    print_progress(lfid, length(all_arrplates),j);
    
end
clear save_mats;

log_printf(lfid, '|\n');

field = 'batchnorm_colsize';


%% Calculate array WT variance
all_arrays = unique(sgadata.arrays);
wild_type_id = strmatch('undefined_sn4757', sgadata.orfnames);

if(isempty(wild_type_id))
    log_printf(lfid, '\n\nTERMINAL WARNING - Cannot calculate array strain variance, no WT screens (%s) found\nWARNING\n', 'undefined_sn4757');
    save('-v7.3',[outputfile,'.mat']);
    return
end

ind2 = query_map{wild_type_id};

array_vars = zeros(length(all_arrays),2);
log_printf(lfid, ['Calculating array WT variance...\n|' blanks(50) '|\n|']);
for i = 1:length(all_arrays)
    ind = intersect(ind2, array_map{all_arrays(i)});
    t = max(sgadata.(field)(ind),1);
    t(isnan(sgadata.(field)(ind))) = NaN;
    array_vars(i,:)=[nanmean(log(t)),nanvar(log(t))];
    
    print_progress(lfid, length(all_arrays),i);
end
log_printf(lfid, '|\n');

nanind = find(isnan(array_vars(:,1)));

% For the NaNs, we need to compute means from the other screens (not just WT).
% This is useful when there's a WT-screen specific issue (e.g., URA3 linkage group).

for i = nanind'
    ind = array_map{all_arrays(i)};
    t = max(sgadata.(field)(ind),1);
    t(isnan(sgadata.(field)(ind))) = NaN;
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
field = 'batchnorm_colsize_nonegs';

log_printf(lfid, ['Computing average for double mutants...\n|' blanks(50) '|\n|']);

for i = 1:length(all_querys)
    
    % BJV this matches Triple WT also
    if all_querys(i) == strmatch('undefined',sgadata.orfnames)
        continue;
    end
    
    ind = query_map{all_querys(i)};
    ind = setdiff(ind, find(isnan(sgadata.(field))));
    
    spots = sgadata.replicateid(ind) + double(sgadata.plateids(ind))*1000;
    [all_spots, m, n] = unique(spots);
    
    % Compute average, variance and number of colonies for each group of spots
    mn = grpstats(log(sgadata.(field)(ind)), spots, 'mean');
    vr = grpstats(log(sgadata.(field)(ind)), spots, 'std').^2;
    nm = grpstats(sgadata.(field)(ind), spots, 'numel');
    
    sgadata.dm_normmean(ind) = mn(n);
    sgadata.dm_normvar(ind) = vr(n);
    sgadata.dm_num(ind) = nm(n);
    
    % Print progress
    print_progress(lfid, length(all_querys), i);
    
end
log_printf(lfid, '|\n');

% Pool across arrayplates for each query
all_arrayplateids = unique(sgadata.arrayplateids);
query_arrplate_vars = [];   %zeros(length(all_querys),length(all_arrplates))+NaN;
query_arrplate_relerr = []; %zeros(length(all_querys),length(all_arrplates))+NaN;

log_printf(lfid, ['Pooling across arrayplates for each query...\n|' blanks(50) '|\n|']);
for i = 1:length(all_querys)
    
    ind = query_map{all_querys(i)};
    
    for j = 1:length(all_arrayplateids)
        
        tmp_ind = find(sgadata.arrayplateids(ind) == all_arrayplateids(j));
        currsets = unique(sgadata.setids(tmp_ind));
        
        for k = 1:length(currsets)
            
            tmp_ind2 = find(sgadata.setids(tmp_ind) == currsets(k));
            curr_ind = tmp_ind(tmp_ind2);
            
            query_arrplate_relerr = [query_arrplate_relerr; i,j,k,sqrt(exp((1./(length(curr_ind)-length(unique(sgadata.arrays(curr_ind)))))*nansum(sgadata.dm_normvar(curr_ind).*((sgadata.dm_num(curr_ind)-1)./sgadata.dm_num(curr_ind))))-1)];
            query_arrplate_vars = [query_arrplate_vars; i,j,k,(1./(length(curr_ind)-length(unique(sgadata.arrays(curr_ind)))))*nansum(sgadata.dm_normvar(curr_ind).*((sgadata.dm_num(curr_ind)-1)./sgadata.dm_num(curr_ind)))];
            
        end
    end
    
    % Print progress
    print_progress(lfid, length(all_querys), i);
end
log_printf(lfid, '|\n');


% Load the single mutant fitness file -----------------------------------------------------------------------------------------------------
%{
[data.f1, data.f2, data.f3] = textread(smfitnessfile,'%s %f %f');
sm_fitness_orfs = data.f1;
sm_fitness = [data.f2 data.f3];
clear data;

[int, a, b] = intersect(sgadata.orfnames, sm_fitness_orfs);
fg_smfit = zeros(length(sgadata.orfnames),2) + NaN;
fg_smfit(a,:) = sm_fitness(b,:);
final_smfit = fg_smfit(:,1);
final_smfit_std = fg_smfit(:,2);
log_printf(lfid, 'Strains not in fitness file: %d\n\n', length(sgadata.orfnames) - length(int));
log_printf(lfid, 'Strains NaN in fitness file: %d\n\n', sum(isnan(fg_smfit(a,1)))); % of those we want (intersection)
%}

% Load the single mutant fitness file -----------------------------------------------------------------------------------------------------
[fitness_data.ORF, fitness_data.SMF, fitness_data.STD] = textread(smfitnessfile,'%s %f %f');
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
    log_printf('\t%s\t: %d\n', fitness_report_header{i}, fitness_report_counts(i));
end
% Load the single mutant fitness file -----------------------------------------------------------------------------------------------------


field = 'batchnorm_colsize';

% replace linkage colony values that got removed before batch correction
if(skip_linkage_mask)
	sgadata.(field)(linkage_cols) = sgadata.filt_jackknife(linkage_cols);
	log_printf(lfid, '* REPLACING linkage colonies.\n');
end

all_arrays = unique(sgadata.arrays);

model_fits = zeros(length(all_arrays),length(sgadata.orfnames)+1) + NaN;
model_fit_std = zeros(length(all_arrays),length(sgadata.orfnames)+1) + NaN;

all_nans = find(isnan(sgadata.(field)));
ind_his3 = strmatch(border_strain_orf, sgadata.orfnames,'exact');

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
    p = polyfit(final_smfit(sgadata.querys(all_ind(ind2))),sgadata.(field)(all_ind(ind2)),1);
    
    curr_data = sgadata.(field)(all_ind);
    
    % Added back (11-03-31)
    curr_data(ind2) = curr_data(ind2) + (nanmean(sgadata.(field)(all_ind(ind2))) - (p(1)*final_smfit(sgadata.querys(all_ind(ind2)))+p(2)));

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

%FIXME why are we short here
sgadata.arraymean_corrected = sgadata.arraymean_corrected';
sgadata.arraymean_corrected(end:length(sgadata.arraymedian)) = nan;
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

all_arrayplateids = unique(sgadata.arrayplateids);
all_arrayplateids_map = zeros(max(all_arrayplateids),1);
all_arrayplateids_map(all_arrayplateids) = 1:length(all_arrayplateids);

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

dm_actual = dm_expected + eps;
dm_actual_std = eps_std;
clear eps eps_std afit;

% Remove "undefined"
if ~exist('skip_wt_remove', 'var')
    skip_wt_remove = false;
end

% remove anything "undefined"
if(~skip_wt_remove)
    ind = strmatch('undefined', sgadata.orfnames);
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
