function [] = define_query_specific_linkage_bjv_wrapper(sga, linkage_filename)
%function [] = define_query_specific_linkage_bjv_wrapper(sga, linkage_filename)
% this is just a wrapper that passes one of my data structures
% to one of Anastasisas tools (which expects her data structure)

sga_ab = struct();
sga_ab.queries = sga.Cannon.Orf(sga.Cannon.isQuery);
sga_ab.arrays  = sga.Cannon.Orf(sga.Cannon.isArray);
sga_ab.scores  = sga.eps(sga.Cannon.isQuery, sga.Cannon.isArray);

[dataset_filt, dataset_corr, dataset_noncorr, lnkg] = define_query_specific_linkage_120121(sga_ab);
cell2csv(linkage_filename, [lnkg.orf num2cell(lnkg.coord_mean)]);

