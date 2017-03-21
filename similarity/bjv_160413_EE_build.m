function bjv_160413_EE_build(sga_exe)

   % load the ts data
   %load ~/Research/SGA/post_review_QC_2016/sga8.mat sga_exe

   % convert struct
   ab_exe = convert_bjv_to_ab(sga_exe);

   % run original code
   cc_ts = ab_code(ab_exe);

   % save the result
   save -v7.3 ~/Research/SGA/sga_living_160429/cc_ts_160429.mat cc_ts

end

function cc_ts_140130 = ab_code(ts_140130)


   %addpath(genpath('Datasets/'))
   %addpath(genpath('Utils/'))
   %javaaddpath('Utils/Java/')

   %load ts_140130_v4

   %% 

   ts_140130.queries_essential = zeros(size(ts_140130.queries));
   ts_140130.arrays_essential = zeros(size(ts_140130.arrays));

   tmp = split_by_delimiter('_', ts_140130.queries);
   inds = find(strncmp('tsq', tmp(:,2),3));
   ts_140130.queries_essential(inds) = 1;

   tmp = split_by_delimiter('_', ts_140130.arrays);
   inds = find(strncmp('tsa', tmp(:,2),3));
   ts_140130.arrays_essential(inds) = 1;

   % Remove DAMPs
   tmp = split_by_delimiter('_', ts_140130.queries);
   inds = find(strncmp('damp', tmp(:,2),4));
   ts_140130.queries(inds) = [];
   ts_140130.queries_gn(inds) = [];
   % ts_140130.queries_orf(inds) = [];
   ts_140130.queries_essential(inds) = [];
   ts_140130.scores_eps(inds,:) = [];
   ts_140130.scores_pvalue(inds,:) = [];
   % ts_140130.cobatch_QQ(inds,:) = [];
   % ts_140130.cobatch_QQ(:,inds) = [];


   % Remove all the y-strains from the dataset
   tmp = split_by_delimiter('_', ts_140130.queries);
   inds = find(strncmp('y', tmp(:,2),1));
   ts_140130.queries(inds) = [];
   ts_140130.queries_gn(inds) = [];
   % ts_140130.queries_orf(inds) = [];
   ts_140130.queries_essential(inds) = [];
   ts_140130.scores_eps(inds,:) = [];
   ts_140130.scores_pvalue(inds,:) = [];
   % ts_140130.cobatch_QQ(inds,:) = [];
   % ts_140130.cobatch_QQ(:,inds) = [];

   %%

   iTsEq = find(ts_140130.queries_essential == 1);
   iTsEa = find(ts_140130.arrays_essential == 1);

   cc1 = MathUtil.computeCorrelation(ts_140130.scores_eps(iTsEq,iTsEa));
   cc1 = cc1+cc1';
   cc2 = MathUtil.computeCorrelation(ts_140130.scores_eps(iTsEq,iTsEa)');
   cc2 = cc2+cc2';


   % Merge QQ and AA
   cc_ts_140130.EE_nodes_gn = unique([ts_140130.queries_gn(iTsEq); ts_140130.arrays_gn(iTsEa)]);
   cc_ts_140130.EE_nodes_Q = zeros(length(cc_ts_140130.EE_nodes_gn),1);
   cc_ts_140130.EE_nodes_A = zeros(length(cc_ts_140130.EE_nodes_gn),1);

   tmp = zeros(length(cc_ts_140130.EE_nodes_gn), length(cc_ts_140130.EE_nodes_gn)) + NaN;
   [t,ind1,ind2] = intersect(cc_ts_140130.EE_nodes_gn, ts_140130.queries_gn(iTsEq));
   tmp(ind1,ind1,1) = cc1(ind2,ind2);
   cc_ts_140130.EE_nodes_Q(ind1) = 1;

   [t,ind1,ind2] = intersect(cc_ts_140130.EE_nodes_gn, ts_140130.arrays_gn(iTsEa));
   tmp(ind1,ind1,2) = cc2(ind2,ind2);
   cc_ts_140130.EE_nodes_A(ind1) = 1;

   tmp(tmp==0) = NaN;

   cc_ts_140130.EE_QQAA = nanmean(tmp,3);
   % cc_ts_140130.EE_nodes_orf = genename2orf(cc_ts_140130.EE_nodes_gn,'orfs');
   cc_ts_140130.EE_nodes_orf = AlleleToOrf(cc_ts_140130.EE_nodes_gn);

   cc_ts_140130.EE_nodes = cc_ts_140130.EE_nodes_gn;   % unique label

end

