function [] = export_final_datasets(fg_merge, ts_merge, dest)
%function [] = export_final_datasets(fg_merge, ts_merge, dest)
   % run this script seperatly on the "living" and "analysis" sets

   % print DAmPs first, then remove them
   file = [dest '/SGA_DAmP.txt'];
   damp = mask_essentials(fg_merge, 'DA');
   print_sga_2015(damp, file);
   damp = mask_essentials(ts_merge, 'DA');
   print_sga_2015(damp, file, false); % append to single file
   fg_merge = mask_essentials(fg_merge, 'dA');
   ts_merge = mask_essentials(ts_merge, 'dA');

   % print the whole dataset
   file = [dest '/SGA_DMA_array.txt'];
   print_sga_2015(fg_merge, file);

   file = [dest '/SGA_TSA_array.txt'];
   print_sga_2015(ts_merge, file);

   % ExE
   file = [dest '/SGA_ExE.txt'];
   [q, a] = E_supp(ts_merge);
   EE = mask_essentials(ts_merge, 'EE');
   EE.Cannon.isQuery(q) = true; % reset essential suppressors with _S tag
   EE.Cannon.isArray(a) = true;
   print_sga_2015(EE, file);

   % NxN
   file = [dest '/SGA_NxN.txt'];
   [q, a] = N_supp(fg_merge);
   NN = mask_essentials(fg_merge, 'NN');
   NN.Cannon.isQuery(q) = true;
   NN.Cannon.isArray(a) = true;
   print_sga_2015(NN, file);

   % NxE ExN
   file = [dest '/SGA_ExN_NxE.txt'];
   [q, ~] = E_supp(fg_merge);
   [~, a] = N_supp(fg_merge);
   EN = mask_essentials(fg_merge, 'EN');
   EN.Cannon.isQuery(q) = true;
   EN.Cannon.isArray(a) = true;
   print_sga_2015(EN, file);

   [q, ~] = N_supp(ts_merge);
   [~, a] = E_supp(ts_merge);
   NE = mask_essentials(ts_merge, 'NE');
   NE.Cannon.isQuery(q) = true;
   NE.Cannon.isArray(a) = true;
   print_sga_2015(NE, file, false);

end


function [isQuery, isArray] = E_supp(sga)
   % return vectors: isQuery & is_essential & is_suppressor
   supp_def_file = '~/SGA/refdata/suppressor_strain_essentiality_160425.csv';
   supp_def = Csv2Cell(supp_def_file);
   E_ix = strcmp('E', supp_def(:,2));
   supp_ess = supp_def(E_ix,1);

   [~, tails] = StripOrfs(sga.Cannon.Orf);
   is_ess = ismember(tails, supp_ess); %is_supp implied
   isQuery = sga.Cannon.isQuery & is_ess;
   isArray = sga.Cannon.isArray & is_ess';
end
function [isQuery, isArray] = N_supp(sga)
   % return vectors: isQuery & is_essential & is_suppressor
   supp_def_file = '~/SGA/refdata/suppressor_strain_essentiality_160425.csv';
   supp_def = Csv2Cell(supp_def_file);
   N_ix = strcmp('N', supp_def(:,2));
   supp_non = supp_def(N_ix,1);

   [~, tails] = StripOrfs(sga.Cannon.Orf);
   is_non = ismember(tails, supp_non); %is_supp implied
   isQuery = sga.Cannon.isQuery & is_non;
   isArray = sga.Cannon.isArray & is_non';
end
