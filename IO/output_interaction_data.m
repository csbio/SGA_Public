%%
% OUTPUT_INTERACTION_DATA - prints out the calculated genetic interactions 
% into a tab delimited text file.
%
%
% Authors: Chad Myers (cmyers@cs.umn.edu), 
% Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
% Benjamin VanderSluis (ben@bjv.us)
%
% Last revision: 2017-12-12
%
%%

function output_interaction_data(outputfile, orfnames, escores, escores_std, ...
   pvals, smfit, smfit_std, dm_expected, dm_actual, dm_actual_std, lfid)

log_printf(lfid, ['Printing output file...\n|' blanks(50) '|\n|']);

fid = fopen([outputfile '.txt'],'w');
% Add column names
% (https://github.com/csbio/SGA_Public/blob/master/Column_Key.md) + 'eps' = Observed_dmf - Expected_dmf
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'Query Orf', 'Array Orf', 'e_score', 'e_score (std)', 'eps score', 'p-value', ... 
    'Query smf', 'Query smf (std)', 'Array smf', 'Array smf (std)', ...
    'Expected dmf', 'Observed dmf', 'Observed dmf (std)');
for i = 1:length(escores)
    for j = 1:length(escores)
        if ~isnan(escores(i,j))
              fprintf(fid,'%s\t%s\t%f\t%f\t%f\t%e\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
                  orfnames{i}, orfnames{j}, escores(i,j), escores_std(i,j), ...
                  dm_actual(i,j) - dm_expected(i,j) ,pvals(i,j), ...
                  smfit(i), smfit_std(i), smfit(j), smfit_std(j), ...
                  dm_expected(i,j), dm_actual(i,j), dm_actual_std(i,j));
        end
    end
    
    print_progress(lfid, length(escores), i);
    
end
fclose(fid);

% output the list of strains, but only if they have scores
% otherwise bad arrays appear
valid_queries = orfnames(sum(~isnan(escores),2)>0);
valid_arrays  = orfnames(sum(~isnan(escores),1)>0);
valid_orfs = union(valid_queries, valid_arrays);
fid = fopen([outputfile '.orf'], 'w');
for i=1:length(valid_orfs)
   fprintf(fid, '%s\n', valid_orfs{i});
end
fclose(fid);

log_printf(lfid,'|\n');
