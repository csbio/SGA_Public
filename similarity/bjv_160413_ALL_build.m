function bjv_160413_ALL_build(fg_merge_reconstruction, ts_merge_reconstruction)
%function bjv_160413_ALL_build(fg_merge_reconstruction, ts_merge_reconstruction)

   % load the ts data
   % load ~/Research/SGA/post_review_QC_2016/merge_reconstructs_nodamp.mat
   % names are _filt instead of _reconstruction if loading from mat

   % convert struct
   ab_nxn = convert_bjv_to_ab(fg_merge_reconstruction);
   ab_exe = convert_bjv_to_ab(ts_merge_reconstruction);

   % run original code
   cc_ALL = ab_code(ab_nxn, ab_exe);

   % save the result
   save -v7.3 ~/Research/SGA/sga_living_160429/cc_ALL_160429.mat cc_ALL



end

function cc_all = ab_code(fg_140130, ts_140130)
   %% Load data

   % addpath(genpath('Datasets/'))
   % addpath(genpath('Utils/'))
   % javaaddpath('Utils/Java/')

   % load fg_140130_v4
   % load ts_140130_v3

   % Remove DAMPs
   tmp = split_by_delimiter('_', fg_140130.queries);
   inds = find(strncmp('damp', tmp(:,2),4));
   fg_140130.queries(inds) = [];
   fg_140130.queries_gn(inds) = [];
   % fg_140130.queries_orf(inds) = [];
   % fg_140130.queries_essential(inds) = [];
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
   % fg_140130.queries_essential(inds) = [];
   fg_140130.scores_eps(inds,:) = [];
   fg_140130.scores_pvalue(inds,:) = [];
   % fg_140130.cobatch_QQ(inds,:) = [];
   % fg_140130.cobatch_QQ(:,inds) = [];

   % Remove DAMPs
   tmp = split_by_delimiter('_', ts_140130.queries);
   inds = find(strncmp('damp', tmp(:,2),4));
   ts_140130.queries(inds) = [];
   ts_140130.queries_gn(inds) = [];
   % ts_140130.queries_orf(inds) = [];
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
   ts_140130.scores_eps(inds,:) = [];
   ts_140130.scores_pvalue(inds,:) = [];
   % ts_140130.cobatch_QQ(inds,:) = [];
   % ts_140130.cobatch_QQ(:,inds) = [];


   %%
   cc_all.nodes_gn = unique([fg_140130.queries_gn; fg_140130.arrays_gn; ts_140130.queries_gn; ts_140130.arrays_gn]);
   % cc_all.nodes_orf = genename2orf(cc_all.nodes_gn,'noannot');
   cc_all.nodes_orf = AlleleToOrf(cc_all.nodes_gn);

   tmp = split_by_delimiter('_', cc_all.nodes_orf);
   cc_all.nodes_orf = tmp(:,1);

   tmp = zeros(length(cc_all.nodes_gn), length(cc_all.nodes_gn),5);

   % Layer 1: FG QQ
   cc_fg.QQ = MathUtil.computeCorrelation(fg_140130.scores_eps);
   cc_fg.QQ = cc_fg.QQ + cc_fg.QQ';

   [~,ind1,ind2] = intersect(fg_140130.queries_gn, cc_all.nodes_gn);
   tmp(ind2,ind2,1) = cc_fg.QQ(ind1,ind1);

   % Layer 2: TS QQ (normalized based on TS AA)
   cc_ts.QQ = MathUtil.computeCorrelation(ts_140130.scores_eps);
   cc_ts.QQ = cc_ts.QQ + cc_ts.QQ';

   [~,iQfgO,iQtsO] = intersect(fg_140130.queries, ts_140130.queries);
   [~,iAfgO,iAtsO] = intersect(fg_140130.arrays, ts_140130.arrays);

   cc1 = MathUtil.computeCorrelation(fg_140130.scores_eps(iQfgO,iAfgO)');  % FG: AA for common arrays across common queries
   cc1 = cc1+cc1';
   cc2 = MathUtil.computeCorrelation(ts_140130.scores_eps(iQtsO,iAtsO)');  % TS: AA for common arrays across common queries
   cc2 = cc2+cc2';

   [data3_norm, table] = table_norm(cc1(:), cc2(:), cc_ts.QQ(:));  % Use the AA to normalize TS QQ
   cc_ts.QQ_norm = NaN(size(cc_ts.QQ));
   cc_ts.QQ_norm(:) = data3_norm;

   [~,ind1,ind2] = intersect(ts_140130.queries_gn, cc_all.nodes_gn);
   tmp(ind2,ind2,2) = cc_ts.QQ_norm(ind1,ind1);

   % Layer 3: FG AA
   cc_fg.AA = MathUtil.computeCorrelation(fg_140130.scores_eps');
   cc_fg.AA = cc_fg.AA + cc_fg.AA';

   [~,ind1,ind2] = intersect(fg_140130.arrays_gn, cc_all.nodes_gn);
   tmp(ind2,ind2,3) = cc_fg.AA(ind1,ind1);

   % Layer 4: TS AA (normalized)

   cc_ts.AA = MathUtil.computeCorrelation(ts_140130.scores_eps');
   cc_ts.AA = cc_ts.AA + cc_ts.AA';

   [data4_norm, table] = table_norm(cc1(:), cc2(:), cc_ts.AA(:)); % Re-use cc1 and cc2 from Layer 2
   cc_ts.AA_norm = NaN(size(cc_ts.AA));
   cc_ts.AA_norm(:) = data4_norm;

   [~,ind1,ind2] = intersect(ts_140130.arrays_gn, cc_all.nodes_gn);
   tmp(ind2,ind2,4) = cc_ts.AA_norm(ind1,ind1);

   % Layer 5: FG-TS AA (normalized)

   iAfgU = setdiff(1:length(fg_140130.arrays), iAfgO);   % FG arrays that are unique
   iAtsU = setdiff(1:length(ts_140130.arrays), iAtsO);   % TS arrays that are unique

   % V1: FG dataset
   fg = fg_140130.scores_eps(iQfgO,:);

   % V2: FG but replace the overlapping array data with data from TS
   fg2 = fg;
   fg2(:,iAfgO) = ts_140130.scores_eps(iQtsO,iAtsO); 

   % V3: union of FG and TS
   fg_ts = [fg_140130.scores_eps(iQfgO,iAfgU) ts_140130.scores_eps(iQtsO, :)];
   fg_ts_arrays = [fg_140130.arrays(iAfgU); ts_140130.arrays];

   cc1 = MathUtil.computeCorrelation(fg');
   cc1 = cc1+cc1';

   cc2 = MathUtil.computeCorrelation(fg2');
   cc2 = cc2+cc2';

   cc3 = MathUtil.computeCorrelation(fg_ts');
   cc3 = cc3+cc3';

   [data3_norm, table] = table_norm(...
       reshape(cc1(iAfgO,iAfgU),[],1), ...
       reshape(cc2(iAfgO,iAfgU),[],1), ...
       reshape(cc3(end-length(ts_140130.arrays)+1:end,1:length(iAfgU)),[],1));

   indr = size(cc3,1)-length(ts_140130.arrays)+1:size(cc3,1);
   indc = 1:length(iAfgU);

   indrm = repmat(indr',1, length(indc));
   indcm = repmat(indc, length(indr),1);

   indi = sub2ind(size(cc3), indrm(:), indcm(:));

   cc3_norm = zeros(size(cc3));
   cc3_norm(indi) = data3_norm;

   cc3_norm = cc3_norm + cc3_norm';

   %fg_ts_arrays_gn = orf2genename(fg_ts_arrays);
   fg_ts_arrays_gn = StrainToAllele(fg_ts_arrays);

   [t,ind1,ind2] = intersect(fg_ts_arrays_gn, cc_all.nodes_gn);
   tmp(ind2,ind2,5) = cc3_norm(ind1,ind1);


   %%
   tmp(tmp==0) = NaN;
   cc_all.QQAA = nanmean(tmp,3);

   % save('Projects/SGA/Cell_map/14-08-27/cc_140130_ALL_140828.mat','cc_all');

end
