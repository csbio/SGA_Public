function [file_smfit, file_smfit_std] = load_smf(smfitnessfile, sgadata, lfid)
% function [file_smfit, file_smfit_std] = load_smf(smfitnessfile, sgadata, lfid)

   smf_fid = fopen(smfitnessfile, 'r');
   fitness_data = struct();
   fitness_data.raw = textscan(smf_fid, '%s%f%f', 'Delimiter', '\t', 'ReturnOnError', false);
   fclose(smf_fid);
   fitness_data.ORF = fitness_data.raw{1};
   fitness_data.SMF = fitness_data.raw{2};
   fitness_data.STD = fitness_data.raw{3};
   fitness_report_header = {'Exact match', 'Partial Match', 'Not Found', 'NaN in file'};
   fitness_report_counts = zeros(1,4);
   fitness_hash = hash_strings(fitness_data.ORF);

   file_smfit = zeros(length(sgadata.orfnames),1) + NaN;
   file_smfit_std = zeros(length(sgadata.orfnames),1) + NaN;
   for i=1:length(sgadata.orfnames)
       if fitness_hash.containsKey(sgadata.orfnames{i})
           file_smfit(i)     = fitness_data.SMF(fitness_hash.get(sgadata.orfnames{i}));
           file_smfit_std(i) = fitness_data.STD(fitness_hash.get(sgadata.orfnames{i}));
           fitness_report_counts(1) = fitness_report_counts(1)+1;
       elseif fitness_hash.containsKey(strip_annotation(sgadata.orfnames{i}, 'last'))
           file_smfit(i)     = fitness_data.SMF(fitness_hash.get(strip_annotation(sgadata.orfnames{i}, 'last')));
           file_smfit_std(i) = fitness_data.STD(fitness_hash.get(strip_annotation(sgadata.orfnames{i}, 'last')));
           fitness_report_counts(2) = fitness_report_counts(2)+1;
       else
           fitness_report_counts(3) = fitness_report_counts(3)+1;
       end
   end
   fitness_report_counts(4) = sum(isnan(file_smfit)) - fitness_report_counts(3);

   log_printf(lfid, 'Fitness file report:\n');
   for i=1:length(fitness_report_header)
       log_printf(lfid, '\t%s\t: %d\n', fitness_report_header{i}, fitness_report_counts(i));
   end
end
