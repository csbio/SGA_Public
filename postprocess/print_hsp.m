function [] = print_hsp(fg_merge_analysis_qvar, ts_merge_analysis_qvar, fg_merge_living, ts_merge_living)

   % restrict "living" data to "analysis" queries
   % did this in post to match Matej's manual effort
   fg_merge_living.Cannon.isArray(~ismember(fg_merge_living.Cannon.Orf, fg_merge_analysis_qvar.Cannon.Orf(fg_merge_analysis_qvar.Cannon.isArray))) = false;
   ts_merge_living.Cannon.isArray(~ismember(ts_merge_living.Cannon.Orf, ts_merge_analysis_qvar.Cannon.Orf(ts_merge_analysis_qvar.Cannon.isArray))) = false;

   F = 'hsp90_data.txt';
   print_hsp_header(F);

   % change the orf labels before printing
   old_targets = {'YEL021W_y13096ura_T26', 'YPL240C_y14083_T26', 'YPL240C_y14084_T26'};
   old_orflabels = {'YEL021W+YDL227C_y13096ura_nat', 'HSP82_y14083', 'YPL240C+YMR186W_y14084'};
   old_commons = {'double mutant control strain', 'hsp82-ts::URA3', 'hsc82^::NAT hsp82-ts::URA3'};

   new_targets = {'YMR186W_sn2174', 'YEL021W_y8835'};
   new_orflabels = {'HSC82_sn2174', 'YEL021W_y8835ura'};
   new_commons = {'hsc82^::NAT', 'single mutant control strain'};

   % old targets from "analysis" data with corrected pvalues
   ix = apply_map(fg_merge_analysis_qvar.Cannon.Map, old_targets);
   fprintf('3:%d\n', length(nonzeros(ix)));
   fg_merge_analysis_qvar.Cannon.isQuery = false(size(fg_merge_analysis_qvar.Cannon.isQuery));
   fg_merge_analysis_qvar.Cannon.isQuery(ix) = true;
   fg_merge_analysis_qvar.Cannon.Orf(ix) = old_orflabels;
   fg_merge_analysis_qvar.Cannon.Common(ix) = old_commons;
   print_sga_append(fg_merge_analysis_qvar, F);

   ix = apply_map(ts_merge_analysis_qvar.Cannon.Map, old_targets);
   fprintf('3:%d\n', length(nonzeros(ix)));
   ts_merge_analysis_qvar.Cannon.isQuery = false(size(ts_merge_analysis_qvar.Cannon.isQuery));
   ts_merge_analysis_qvar.Cannon.isQuery(ix) = true;
   ts_merge_analysis_qvar.Cannon.Orf(ix) = old_orflabels;
   ts_merge_analysis_qvar.Cannon.Common(ix) = old_commons;
   print_sga_append(ts_merge_analysis_qvar, F);

   % new targets from "living" data 
   ix = apply_map(fg_merge_living.Cannon.Map, new_targets);
   fprintf('2:%d\n', length(nonzeros(ix)));
   fg_merge_living.Cannon.isQuery = false(size(fg_merge_living.Cannon.isQuery));
   fg_merge_living.Cannon.isQuery(ix) = true;
   fg_merge_living.Cannon.Orf(ix) = new_orflabels;
   fg_merge_living.Cannon.Common(ix) = new_commons;
   print_sga_append(fg_merge_living, F);

   ix = apply_map(ts_merge_living.Cannon.Map, new_targets);
   fprintf('2:%d\n', length(nonzeros(ix)));
   ts_merge_living.Cannon.isQuery = false(size(ts_merge_living.Cannon.isQuery));
   ts_merge_living.Cannon.isQuery(ix) = true;
   ts_merge_living.Cannon.Orf(ix) = new_orflabels;
   ts_merge_living.Cannon.Common(ix) = new_commons;
   print_sga_append(ts_merge_living, F);
end

function [] = print_hsp_header(filename)
   fid = fopen(filename, 'w');
   header = {
   'Query Strain ID',
   'Query description',
   'Array Strain ID',
   'Array description',
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
   fclose(fid);
end

function [] = print_sga_append(sga, filename)
   % replace array common names with alleles
   fid = fopen('~/Research/Data/YeastGeneMap/strain_orf_common_allele_151109.txt', 'r');
   A = textscan(fid, '%s%s%s%s', 'Delimiter', '\t', 'ReturnOnError', false);
   fclose(fid);

   map = hash_strings(A{1});
   array_ids = apply_map(map, sga.Cannon.Orf(sga.Cannon.isArray));

   % start with common to catch misses
   array_alleles = sga.Cannon.Common(sga.Cannon.isArray);
   % replace with alleles
   array_alleles(array_ids ~=0) = A{4}(nonzeros(array_ids));
   % put into place
   sga.Cannon.Common(sga.Cannon.isArray) = array_alleles;

   fid = fopen(filename, 'a');
   % print in row-major order
   % for/find loops only work on row-vectors
   % so let's check assumptions
   assert(size(sga.Cannon.isQuery, 2) == 1);
   assert(size(sga.Cannon.isArray, 1) == 1);

   q_ix = find(sga.Cannon.isQuery');
   a_ix = find(sga.Cannon.isArray);

   for q=q_ix
      for a=a_ix
         if isnan(sga.eps(q,a))
            continue;
         end

         % first 5 text columns
         fprintf(fid, '%s\t%s\t%s\t%s\t%s\t', ...
            sga.Cannon.Orf{q}, ...
            sga.Cannon.Common{q}, ...
            sga.Cannon.Orf{a}, ...
            sga.Cannon.Common{a}, ...
            sga.source_labels{sga.source(q,a)});
         % score
         fprintf(fid, '%.3f\t%.3e\t', sga.eps(q,a), sga.pvl(q,a));                                                

         % query & array fitnesses
         fprintf(fid, '%.3f\t%.3f\t', sga.qfit(q,a), sga.afit(q,a));                                                                        

         % double
         fprintf(fid, '%.3f\t%.3f\n', sga.dbl(q,a), sga.dbl_std(q,a)); 
      end
   end
   fclose(fid);
end
