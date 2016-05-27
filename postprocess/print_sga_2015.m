function[] = print_sga_2015(sga, filename, newfile)
   if ~exist('newfile', 'var')
      newfile=true;
   end

   alleles = StrainToAllele(sga.Cannon.Orf);
   if(newfile)
      fid = fopen(filename, 'w');
      header = {
      'Query Strain ID', 
      'Query allele name', 
      'Array Strain ID', 
      'Array allele name', 
      'Arraytype/Temp', 
      'Genetic interaction score (ε)', 
      'P-value', 
      'Query single mutant fitness (SMF)', 
      'Array SMF', 
      'Double mutant fitness', 
      'Double mutant fitness standard deviation'};

      h_line = header{1};
      for i=2:length(header)
         h_line = [h_line sprintf('\t%s', header{i})];
      end
      fprintf(fid, '%s\n', h_line);
   else
      fid = fopen(filename, 'a');
   end

   % print in row-major order
   % for/find loops only work on row-vectors
   % so let's check assumptions
   assert(size(sga.Cannon.isQuery, 2) == 1);
   assert(size(sga.Cannon.isArray, 1) == 1);

   q_ix = find(sga.Cannon.isQuery');
   a_ix = find(sga.Cannon.isArray);


   for q=q_ix
      for a=a_ix
         % skip NaNs
         if isnan(sga.eps(q,a))
            continue;
         end

         % first 5 text columns
         fprintf(fid, '%s\t%s\t%s\t%s\t%s\t', ...
            sga.Cannon.Orf{q}, ...
            alleles{q}, ...
            sga.Cannon.Orf{a}, ...
            alleles{a}, ...
            sga.source_labels{sga.source(q,a)});

         % score
         fprintf(fid, '%.4f\t%.3e\t', sga.eps(q,a), sga.pvl(q,a));

         % query & array fitnesses
         fprintf(fid, '%.4f\t%.4f\t', ...
            sga.qfit(q,a), sga.afit(q,a));

         % double
         fprintf(fid, '%.4f\t%.4f\n', sga.dbl(q,a), sga.dbl_std(q,a));
      end
   end
   fclose(fid);
end
