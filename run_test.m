
addpath('~/SGA/Main/corrections');
addpath('~/SGA/Main/IO');
addpath('~/SGA/Main/util');
linkage_file = 'refdata/linkage_estimate_curated_151111.txt';
coord_file = 'refdata/chrom_coordinates_150617.tab';
lfid = -11;


%load scored_sga_ts_t30_150617_scored_150904.mat sgadata all_querys all_arrays query_map array_map wild_type
%[all_linkage_cols, non_spec] = filter_linkage_colonies(sgadata, linkage_file, ...
%                               coord_file, all_querys, all_arrays, query_map, array_map, ...
%                               wild_type, lfid);



% load each of the four latest datasets in turn, and print out linkage receipts
base = '/home/grad06/bvander/SGA/Main/scored/151015/';
mats = {'scored_sga_fg_t26_151015_scored_151030.mat',...
         'scored_sga_fg_t30_151015_scored_151030.mat',...
         'scored_sga_ts_t26_151015_scored_151030.mat',...
         'scored_sga_ts_t30_151015_scored_151030.mat'};

vars = 'sgadata all_querys all_arrays query_map array_map wild_type';

for i=1:length(mats)
   eval(sprintf('load %s%s %s', base, mats{i}, vars));
   [all_linkage_cols, non_spec] = filter_linkage_colonies(sgadata, linkage_file, ...
                                  coord_file, all_querys, all_arrays, query_map, array_map, ...
                                  wild_type, lfid);
end



