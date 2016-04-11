function [sga] = add_fitnesses(sga)

   file_26 = 'refdata/smf_t26_130417.txt';
   fid = fopen(file_26, 'r');
   A = textscan(fid, '%s%f%f', 'Delimiter', '\t', 'ReturnOnError', false);
   fclose(fid);
   hash_26 = java.util.HashMap(length(A{1}));
   for i=1:length(A{1})
      hash_26.put(java.lang.String(A{1}{i}), A{2}(i));
   end
   fit26 = NaN(sga.Cannon.GENES,1);
   for i=1:sga.Cannon.GENES
      f = hash_26.get(sga.Cannon.Orf{i});
      if ~isempty(f)
         fit26(i)=f;
      end
   end

   file_30 = 'refdata/smf_t30_130417.txt';
   fid = fopen(file_30, 'r');
   A = textscan(fid, '%s%f%f', 'Delimiter', '\t', 'ReturnOnError', false);
   fclose(fid);
   hash_30 = java.util.HashMap(length(A{1}));
   for i=1:length(A{1})
      hash_30.put(java.lang.String(A{1}{i}), A{2}(i));
   end
   fit30 = NaN(sga.Cannon.GENES,1);
   for i=1:sga.Cannon.GENES
      f = hash_30.get(sga.Cannon.Orf{i});
      if ~isempty(f)
         fit30(i)=f;
      end
   end


   qfit_30 = repmat(fit30, [1,size(sga.eps,2)]);
   qfit_26 = repmat(fit26, [1,size(sga.eps,2)]);

   afit_30 = repmat(fit30', [size(sga.eps,1),1]);
   afit_26 = repmat(fit26', [size(sga.eps,1),1]);


   sga.qfit = NaN(size(sga.eps));
   sga.afit = NaN(size(sga.eps));

   mask_26 = (sga.source == 1) | (sga.source == 3);
   mask_30 = (sga.source == 2) | (sga.source == 4);


   sga.qfit(mask_26) = qfit_26(mask_26);
   sga.qfit(mask_30) = qfit_30(mask_30);
   sga.afit(mask_26) = afit_26(mask_26);
   sga.afit(mask_30) = afit_30(mask_30);


end




