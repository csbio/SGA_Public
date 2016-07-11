% test_30
% Check individual interactions

% load the txt files, they will be left in the workspace and saved in a matfile
% in case further testing is required.
sga_std = load_sga_epsilon_from_scorefile('Sci_Batch40_standard.txt', 'Sci_Batch40_standard.orf');
sga_out = load_sga_epsilon_from_scorefile('Sci_Batch40_output.txt', 'Sci_Batch40_output.orf');
save result_30.mat sga_std sga_out

% check the total number of valid scores
eps = sga_std.eps(sga_std.Cannon.isQuery, sga_std.Cannon.isArray);
pvl = sga_std.pvl(sga_std.Cannon.isQuery, sga_std.Cannon.isArray);
valid_std = sum(sum(~isnan(eps) & ~isnan(pvl)));

eps = sga_out.eps(sga_out.Cannon.isQuery, sga_out.Cannon.isArray);
pvl = sga_out.pvl(sga_out.Cannon.isQuery, sga_out.Cannon.isArray);
valid_out = sum(sum(~isnan(eps) & ~isnan(pvl)));

if valid_std == valid_out
   fprintf('Number of scores: Match. (%d)\n', valid_std);
elseif valid_out > valid_std
   fprintf('Number of scores: INCREASED (%d + %d)\n', valid_std, valid_out-valid_std);
   fprintf('Missing values will not count further as mismatches.\n');
else
   fprintf('Number of scores: DECREASED (%d - %d)\n', valid_std, valid_std-valid_out);
   fprintf('Missing values will not count further as mismatches.\n');
end
clear eps pvl valid_std valid_out 

% intersect the data, and examine scores themselves
[E_std, E_out, P_std, P_out] = eps_intersect(sga_std, sga_out);
E_std = E_std(:);
E_out = E_out(:);
P_std = P_std(:);
P_out = P_out(:);
valid = ~isnan(E_std) & ~isnan(E_out) & ~isnan(P_std) & ~isnan(P_out);

% Epsilon: mismatch
miss = E_std(valid) ~= E_out(valid);
fprintf('Epsilon mismatches: %d / %d (%d%%)\n',...
       sum(miss), length(miss), round(100*sum(miss)/length(miss)));

% Epsilon: error
abs_error = abs(E_std(valid) - E_out(valid));
fprintf('Epsilon mean absolute difference: %.2e\n', mean(abs_error));
fprintf('Epsilon std  absolute difference: %.2e\n', std(abs_error));

% Epsilon correlation 
[r, p] = corr([E_std(valid) E_out(valid)], 'rows', 'pairwise', 'type', 'Pearson');
fprintf('Epsilon Pearson correlation R: %.2f  P: %.2e\n', r(1,2), p(1,2));

% Pvalue: mismatch
miss = P_std(valid) ~= P_out(valid);
fprintf('Pvalue mismatches: %d / %d (%d%%)\n',...
       sum(miss), length(miss), round(100*sum(miss)/length(miss)));

% Epsilon / pvalue scatter plots
% limit scatter plot to 10K sampleing
valid = random_subset(valid, 10000);
subplot(1,2,1)
scatter(E_std(valid), E_out(valid), 5, '.');
xlabel('Epsilon Standard');
ylabel('Epsilon Output');

subplot(1,2,2)
scatter(P_std(valid), P_out(valid), 5, '.');
xlabel('Pvalue Standard');
ylabel('Pvalue Output');

saveas(gcf(), 'result_30.pdf');
clear E_std E_out P_std P_out miss valid abs_error r p

