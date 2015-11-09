function [] = export_merged_temperature(sga_outputfile)
%function [] = export_merged_temperature(sga_outputfile)
%
% merges temperatures, saves results, clusters results
% call me just once, with any of the four outputfiles

   dirname = split_by_delimiter('/', sga_outputfile);                                                                                                                                
   basename= split_by_delimiter('_', dirname{end});

   int_dirname = [join_by_delimiter(dirname(1:end-1), '/') '/interactions/'];
   clus_dirname = [join_by_delimiter(dirname(1:end-1), '/') '/clustergrams/']; 

   eval(sprintf('load %s/sga_fg_t30.mat', int_dirname));
   eval(sprintf('load %s/sga_fg_t26.mat', int_dirname));
   eval(sprintf('load %s/sga_ts_t30.mat', int_dirname));
   eval(sprintf('load %s/sga_ts_t26.mat', int_dirname));

   fg_merge = merge_temperatures(sga_fg_t30, 'FG30', sga_fg_t26, 'FG26');
   eval(sprintf('save -v7.3 %s/fg_merge.mat fg_merge', int_dirname));

   ts_merge = merge_temperatures(sga_ts_t26, 'TS26', sga_ts_t30, 'TS30');
   eval(sprintf('save -v7.3 %s/ts_merge.mat ts_merge', int_dirname));

   basename{1} = 'clustergram'; 
   basename{4} = 'merge';

   basename{3} = 'fg';
   clus_basename = join_by_delimiter(basename, '_');                                                                                                                            
   generate_fg_clustergram(fg_merge, [clus_dirname clus_basename]); 

   basename{3} = 'ts';
   clus_basename = join_by_delimiter(basename, '_');                                                                                                                            
   generate_fg_clustergram(ts_merge, [clus_dirname clus_basename]); 

end
