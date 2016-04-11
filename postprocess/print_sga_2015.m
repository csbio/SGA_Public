function[] = print_sga_2015(sga, filename, newfile)
   if ~exist('newfile', 'var')
      newfile=true;
   end


   % prep an allele mapping
   allele_file = '~/Research/Data/YeastGeneMap/strain_orf_common_allele_151109.txt';
   fid = fopen(allele_file, 'r');
   A = textscan(fid, '%s%s%s%s', 'Delimiter', '\t', 'ReturnOnError', false);
   fclose(fid);
   a_map = java.util.HashMap();
   for i=1:length(A{1})
      a_map.put(java.lang.String(A{1}{i}), java.lang.String(A{4}{i}));
   end

   alleles = sga.Cannon.Common;
   for i=1:length(alleles)
      a = a_map.get(sga.Cannon.Orf{i});
      if ~isempty(a)
         alleles{i} = a;
      end
   end
   
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
