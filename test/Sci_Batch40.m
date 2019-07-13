% 40 randomly selected batches from the 2010 dataset
% (supplementary data file 4)
% With strain IDs added later by BJV

% Sci_40.txt.gz summary:
% type    sn    tsq    damp    dma    tsa    y    unann    total
% query   122   8      7       0      0      0    0        137
% array   0     0      0       4293   0      0    0        4293

% note that this value for wild_type is only valid for old datasets, 
% the value used in new data is 'URA3control_sn4757'
wild_type = 'undefined_sn4757'; 

% I can vouch for this number... Totally random.
random_seed = 42;

skip_perl_step = false;
border_strain_orf = 'YOR202W_dma1';
smfitnessfile = 'refdata/smf_t30_150803.txt';
linkagefile = 'refdata/linkage_estimate_curated_160426.txt';
removearraylist = 'refdata/bad_strains_160303.csv';
coord_file = 'refdata/chrom_coordinates_150617.tab';

inputfile = 'test/Sci_Batch40/Sci_Batch40.txt';
outputfile = 'test/Sci_Batch40/Sci_Batch40_output';

skip_batch_correction = true;

compute_sgascore

