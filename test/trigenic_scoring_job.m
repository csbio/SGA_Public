%% Step 1: Parameter Setup

% A. Mandatory Parameters: provides locations of requried files
inputfile = 'test/trigenic_scoring_test/raw_triple_test_small.txt'; % unzip raw_triple_test_small.txt.gz
outputfile= 'test/trigenic_scoring_test/scored_triple_test_small_digenic.txt';

smfitnessfile = 'refdata/smf_DmQueryStd_MiniArray_180208.tsv'; 
linkagefile = 'refdata/linkage_estimate_curated_160426.txt'; 
coord_file = 'refdata/chrom_coordinates_150617.tab'; 
removearraylist = 'refdata/bad_strains_160303.csv'; 

% B. Overridable parameters
border_strain_orf = 'YOR202W_dma1';
wild_type = 'URA3control+YDL227C_y13096';
skip_perl_step = false;

% C. Optional parameters
% skip_linkage_detection = false;
% skip_linkage_mask = false;
% skip_wt_remove = false;
random_seed = 42; % Reproducibility


%% Step 2: Add the SGA repository into path

% Assuming the current directory contains base_dir (relative path); if not, 
% use the absolute path of the base_dir here
base_dir = 'SGA_Public-master'; 
addpath(base_dir);

% Add necessary subdirectories
add_SGAPATH();


%% Step 3: Generate SGA scores
cd(get_SGAROOT); % Make sure you are inside the base_dir
compute_sgascore


%% Step 4-7: Generate Tau-SGA scores

% Step 4: Load from standard sga output file (12 columns)
sga = load_sga_epsilon_from_scorefile([outputfile '.txt'], [outputfile '.orf']);

% Step 5: Provide an assignent file (inside SGA_Public-master/refdata)
assignments = 'assignment_file_170328.csv'; 

% Step 6: Calculate Tau scores (trigenic scores)
sga_triple = score_trigenic_interactions(sga, assignments);

% Step 7: Print/write the Tau scores
tau_output_file = strrep(outputfile, 'digenic', 'trigenic');
print_trigenic(sga, sga, tau_output_file);


%% Step 8: Clustering
generate_fg_clustergram(sga_triple, tau_output_file);
