% 40 randomly selected batches from the 2010 dataset
% (supplementary data file 4)
% With strain IDs added later by BJV

% Sci_40.txt.gz summary:
% type    sn    tsq    damp    dma    tsa    y    unann    total
% query   122   8      7       0      0      0    0        137
% array   0     0      0       4293   0      0    0        4293

% note that this value for wild_type is only valid for old datasets, 
% new value is 'URA3control_sn4757'
wild_type = 'undefined_sn4757'; 

skip_perl_step = false;
border_strain_orf = 'YOR202W_dma1';
smfitnessfile = 'refdata/smf_t30_150803.txt';
linkagefile = 'refdata/linkage_estimate_curated_160426.txt';
removearraylist = 'refdata/bad_strains_160303.csv';
coord_file = 'refdata/chrom_coordinates_150617.tab';

inputfile = 'rawdata/Sci_Batch40/Sci_Batch40.txt';
outputfile = 'scored/Sci_Batch40/Sci_Batch40';

compute_sgascore

