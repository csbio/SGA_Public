function bjv_160413_NN_build(sga_nxn) 
%function bjv_160413_NN_build(sga_nxn) 

   % load the ts data
   % load ~/Research/SGA/post_review_QC_2016/sga8.mat sga_nxn

   % convert struct
   ab_nxn = convert_bjv_to_ab(sga_nxn);

   % run original code
   cc_fg = ab_code(ab_nxn);

   % save the result
   save -v7.3 ~/Research/SGA/sga_living_160429/cc_fg_160429.mat cc_fg

end

function cc_fg_140130 = ab_code(fg_140130)
   % addpath(genpath('Datasets/'))
   % addpath(genpath('Utils/'))
   % javaaddpath('Utils/Java/')

   % load fg_140130_v4

   %%

   fg_140130.queries_essential = zeros(size(fg_140130.queries));
   fg_140130.arrays_essential = zeros(size(fg_140130.arrays));

   tmp = split_by_delimiter('_', fg_140130.queries);
   inds = find(strncmp('tsq', tmp(:,2),3));
   fg_140130.queries_essential(inds) = 1;

   % Remove DAMPs
   tmp = split_by_delimiter('_', fg_140130.queries);
   inds = find(strncmp('damp', tmp(:,2),4));
   fg_140130.queries(inds) = [];
   fg_140130.queries_gn(inds) = [];
   % fg_140130.queries_orf(inds) = [];
   fg_140130.queries_essential(inds) = [];
   fg_140130.scores_eps(inds,:) = [];
   fg_140130.scores_pvalue(inds,:) = [];
   % fg_140130.cobatch_QQ(inds,:) = [];
   % fg_140130.cobatch_QQ(:,inds) = [];

   % Remove all the y-strains from the dataset
   tmp = split_by_delimiter('_', fg_140130.queries);
   inds = find(strncmp('y', tmp(:,2),1));
   fg_140130.queries(inds) = [];
   fg_140130.queries_gn(inds) = [];
   % fg_140130.queries_orf(inds) = [];
   fg_140130.queries_essential(inds) = [];
   fg_140130.scores_eps(inds,:) = [];
   fg_140130.scores_pvalue(inds,:) = [];
   % fg_140130.cobatch_QQ(inds,:) = [];
   % fg_140130.cobatch_QQ(:,inds) = [];


   %%
   iFgNq = find(fg_140130.queries_essential == 0);
   iFgNa = find(fg_140130.arrays_essential == 0);

   cc1 = MathUtil.computeCorrelation(fg_140130.scores_eps(iFgNq,iFgNa));
   cc1 = cc1+cc1';
   cc2 = MathUtil.computeCorrelation(fg_140130.scores_eps(iFgNq,iFgNa)');
   cc2 = cc2+cc2';


   % Merge QQ and AA
   % AB wrote this for the natural NN ORF space, I'm converting to allele space for consistancy
   % with the EE code
   % cc_fg_140130.NN_nodes = unique([fg_140130.queries_orf(iFgNq); fg_140130.arrays_orf(iFgNa)]);
   cc_fg_140130.NN_nodes_gn = unique([fg_140130.queries_gn(iFgNq); fg_140130.arrays_gn(iFgNa)]);
   cc_fg_140130.NN_nodes_Q = zeros(length(cc_fg_140130.NN_nodes_gn),1);
   cc_fg_140130.NN_nodes_A = zeros(length(cc_fg_140130.NN_nodes_gn),1);

   tmp = zeros(length(cc_fg_140130.NN_nodes_gn), length(cc_fg_140130.NN_nodes_gn)) + NaN;
   % [t,ind1,ind2] = intersect(cc_fg_140130.NN_nodes, fg_140130.queries_orf(iFgNq));
   [t,ind1,ind2] = intersect(cc_fg_140130.NN_nodes_gn, fg_140130.queries_gn(iFgNq));
   tmp(ind1,ind1,1) = cc1(ind2,ind2);
   cc_fg_140130.NN_nodes_Q(ind1) = 1;

   % [t,ind1,ind2] = intersect(cc_fg_140130.NN_nodes, fg_140130.arrays_orf(iFgNa));
   [t,ind1,ind2] = intersect(cc_fg_140130.NN_nodes_gn, fg_140130.arrays_gn(iFgNa));
   tmp(ind1,ind1,2) = cc2(ind2,ind2);
   cc_fg_140130.NN_nodes_A(ind1) = 1;

   tmp(tmp==0) = NaN;

   cc_fg_140130.NN_QQAA = nanmean(tmp,3);

   % cc_fg_140130.NN_nodes_gn = orf2genename(cc_fg_140130.NN_nodes);
   cc_fg_140130.NN_nodes_orf = AlleleToOrf(cc_fg_140130.NN_nodes_gn);
   cc_fg_140130.NN_nodes = cc_fg_140130.NN_nodes_gn;

   % Quick fix: strip the wrong annotations accidentally appended to the genenames
   % as I'm not using AB's name mapper, I should not inherit any bugs
   % tmp = split_by_delimiter('_', cc_fg_140130.NN_nodes_gn);
   % cc_fg_140130.NN_nodes_gn = tmp(:,1);
   % tmp = split_by_delimiter('-', cc_fg_140130.NN_nodes_gn);
   % inds = find(~cellfun(@isempty, tmp(:,2)) & ~strcmp(upper(cc_fg_140130.NN_nodes_gn), cc_fg_140130.NN_nodes));
   % cc_fg_140130.NN_nodes_gn(inds) = tmp(inds,1);

   % save('Projects/SGA/Cell_map/14-08-27/cc_fg_140130_NN_140827.mat','cc_fg_140130');



end
