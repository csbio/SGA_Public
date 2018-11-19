function[] = print_tau_2018(sga_raw, sga_adj, filename, int_only)
%function[] = print_tau_2018(sga_raw, sga_adj, filename, int_only)
% if int_only: print only significant digenics (+/-) and sig trigen (- only)

   fprintf('skipping select_set\n');
   % sga_adj = trigenic_select_set(sga_adj, {'PPI controls', 'Tspace'});
   sga_raw.Cannon.isQuery = sga_adj.Cannon.isQuery;
   sga_raw.Cannon.isArray = sga_adj.Cannon.isArray;

   if int_only
      % annotate all interactions as digenic, novel trigenic, 
      % selection doesn't matter, just need the map
      [map] = trigenic_select_interactions6(sga_adj);
   end

   alleles = trigenic_strain_to_allele(sga_adj.Cannon.Orf);
   fid = fopen(filename, 'w');
   header = {
   'Query strain ID',  %
   'Query allele name',  %
   'Array strain ID',  %
   'Array allele name',  %
   'Combined mutant type',  %
   'Raw genetic interaction score (epsilon)', 
   'Adjusted genetic interaction score (epsilon or tau)' %
   'P-value',  %
   'Query single/double mutant fitness', 
   'Array single mutant fitness', 
   'Double/triple mutant fitness', 
   'Double/triple mutant fitness standard deviation'};

   if int_only
      header = [header([1:5, 7:8]); {'Digenic, Modified trigenic, or Novel trigenic'}];
   end

   h_line = header{1};
   for i=2:length(header)
      h_line = [h_line sprintf('\t%s', header{i})];
   end
   fprintf(fid, '%s\n', h_line);

   % print in row-major order
   % for/find loops only work on row-vectors
   % so let's check assumptions
   assert(size(sga_adj.Cannon.isQuery, 2) == 1);
   assert(size(sga_adj.Cannon.isArray, 1) == 1);

   q_ix = find(sga_adj.Cannon.isQuery');
   a_ix = find(sga_adj.Cannon.isArray);

   [qshrtL,qshrtR] = SplitOrfs(strip_annotation(sga_adj.Cannon.Orf));
   HO_map = strcmp('YDL227C', [qshrtL, qshrtR]);

   for q=q_ix
      for a=a_ix
         % skip NaNs
         if isnan(sga_adj.eps(q,a)) || isnan(sga_adj.pvl(q,a))
            continue;
         end
         
         % skip non-interactions if int only
         if int_only
            % digencis
            if sum(HO_map(q,:)) > 0
               if ((sga_adj.eps(q,a) > -0.08) && (sga_adj.eps(q,a) < 0.08)) || (sga_adj.pvl(q,a) >= 0.05)
                  continue
               end
            % trigenics
            else
               % if (sga_adj.eps(q,a) > -0.08) || (sga_adj.pvl(q,a) >= 0.05)
               if ((sga_adj.eps(q,a) > -0.08) && (sga_adj.eps(q,a) < 0.08)) || (sga_adj.pvl(q,a) >= 0.05)
                  continue
               end
            end
         end

         % skip trigenic scores > 0 in either file
         % trigenics
         % if sum(HO_map(q,:)) == 0
            % if (sga_adj.eps(q,a) > 0)
               % continue
            % end
         % end

         % first 4 text columns
         fprintf(fid, '%s\t%s\t%s\t%s\t', ...
            sga_adj.Cannon.Orf{q}, ...
            alleles{q}, ...
            sga_adj.Cannon.Orf{a}, ...
            alleles{a});

         % query type
         if sum(HO_map(q,:)) > 0
            fprintf(fid, 'digenic\t');
         else
            fprintf(fid, 'trigenic\t');
         end

         % scores
         if int_only 
            fprintf(fid, '%.6f\t%.3e\t', ...
                 sga_adj.eps(q,a), sga_adj.pvl(q,a));
            fprintf(fid, '%s\n', map{q,a});
         else
            fprintf(fid, '%.6f\t%.6f\t%.3e\t', sga_raw.eps(q,a), ...
                 sga_adj.eps(q,a), sga_adj.pvl(q,a));

            % query & array fitnesses
            fprintf(fid, '%.4f\t%.4f\t', ...
               sga_adj.qfit(q,a), sga_adj.afit(q,a));

            % "double/trip fitness and std-dev
            fprintf(fid, '%.4f\t%.4f\n', ...
            sga_adj.dbl(q,a), sga_adj.dbl_std(q,a));
         end

      end
   end

   fclose(fid);
end
