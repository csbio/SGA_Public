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
   EE = mask_essentials(ts_merge, 'EE');
   print_sga_2015(EE, file);

   % NxN
   file = [dest '/SGA_NxN.txt'];
   NN = mask_essentials(fg_merge, 'NN');
   print_sga_2015(NN, file);

   % NxE ExN
   file = [dest '/SGA_ExN_NxE.txt'];
   EN = mask_essentials(fg_merge, 'EN');
   print_sga_2015(EN, file);

   NE = mask_essentials(ts_merge, 'NE');
   print_sga_2015(NE, file, false);

end
