function [] = compare_sga_structs_recall(sga, labels, standard, std_label)
%function[] = compare_sga_structs(_recallsga, labels, standard, std_label)

	PLOTS = true;
	SINGLE_PLOT = false;

	if(PLOTS)
		if(SINGLE_PLOT)
			neg_handle = subplot(1,2,1);
			pos_handle = subplot(1,2,2);
		else
			neg_handle = figure();
			pos_handle = figure();
		end
	end
	CM = colormap();

	
	for i=1:length(sga)
		%remove nans (not found in std)
		% mask out any queries or arrays not in the standard
		sga{i}.Cannon.isQuery(~ismember(strip_annotation(sga{i}.Cannon.Orf), standard.orfs)) = false;
		sga{i}.Cannon.isArray(~ismember(strip_annotation(sga{i}.Cannon.Orf), standard.orfs)) = false;

		% perform the analysis
      	PR_neg = help_curve_the_pr(sga{i}, standard, true);
	    PR_pos = help_curve_the_pr(sga{i}, standard, false);

		% draw the results
		if(PLOTS)
			% background, draw this for each SGA in case they are different.
			% background is the same for pos and negative.
			colors{i} = CM(floor(size(CM,1)/(length(sga)+1))*i,:);
			[background_rate(i), background_count(i)] = compute_background_rate(standard, sga{i});

			% determine and print percent recall
			fprintf('%s recall rate (neg): %.2f\n', labels{i}, full(PR_neg(end,3)) / full(background_count(i)));
			fprintf('%s recall rate (pos): %.2f\n', labels{i}, full(PR_pos(end,3)) / full(background_count(i)));

			if(SINGLE_PLOT)
				subplot(neg_handle);
			else
				figure(neg_handle);
			end
			semilogx(PR_neg(:,3), PR_neg(:,2), 'Color', colors{i}, 'LineWidth', 2);
			xx = get(gca, 'XLim');
			set(gca, 'XLim', [min(xx(2)-1, 10), xx(2)]);
			hold on;

			if(SINGLE_PLOT)
				subplot(pos_handle);
			else
				figure(pos_handle);
			end

			semilogx(PR_pos(:,3), PR_pos(:,2), 'Color', colors{i}, 'LineWidth', 2);
			xx = get(gca, 'XLim');
			set(gca, 'XLim', [min(xx(2)-1, 10), xx(2)]);
			hold on;
		end
	end

	if(PLOTS)
		if(SINGLE_PLOT)
			subplot(neg_handle);
		else
			figure(neg_handle);
		end
		legend(labels, 'Interpreter', 'None');
		% Plot the background_rate afterward so they don't show up in the legend
		for i=1:length(sga)
			xx = get(gca, 'XLim');
			plot(xx, [background_rate(i) background_rate(i)], ...
				'LineWidth', 2, 'LineStyle', '--', 'Color', colors{i});
		end
		xlabel('Recall');
		ylabel('Precision');
		title(['Negatives - ' std_label], 'Interpreter', 'None');
		hold off

		if(SINGLE_PLOT)
			subplot(pos_handle);
		else
			figure(pos_handle);
		end
		legend(labels, 'Interpreter', 'None');
		for i=1:length(sga)
			xx = get(gca, 'XLim');
			plot(xx, [background_rate(i) background_rate(i)], ...
				'LineWidth', 2, 'LineStyle', '--', 'Color', colors{i});
		end
		xlabel('Recall');
		ylabel('Precision');
		title(['Positives - ' std_label], 'Interpreter', 'None');
		hold off
	end
end

function [PR] = help_curve_the_pr(sga, standard, negativesTF)
	EPS = sga.eps;
	EPS(sga.pvl > 0.05) = 0;

	if(negativesTF)
   		EPS(EPS > -0.08) = 0;
	else
   		EPS(EPS < 0.08) = 0;
	end
	EPS(isnan(EPS)) = 0;

	%optional exclude stringent
	%EPS(abs(EPS)>0.12) = 0;

   EPS = EPS(sga.Cannon.isQuery, sga.Cannon.isArray);
   
   [r, c, v] = find(EPS);
   Qname = strip_annotation(sga.Cannon.Orf(sga.Cannon.isQuery), 'first');
   Aname = strip_annotation(sga.Cannon.Orf(sga.Cannon.isArray), 'first');
   
   [tmp, ix] = sort(abs(v), 'descend');
   r = r(ix);
   c = c(ix);
   v = v(ix);

   [PR, ix] = help_eval_ppi_pairs(standard, [Qname(r) Aname(c)]);
   PR = [v(ix) PR];



end
function [PR, ix] = help_eval_ppi_pairs(standard, pair_list)

	name_map = hash_strings(standard.orfs);
   	rows = apply_map(name_map, pair_list(:,1));
   	cols = apply_map(name_map, pair_list(:,2));

	linear_ind = sub2ind(size(standard.matrix), rows, cols);
	linear_vals = standard.matrix(linear_ind);
	if(islogical(standard.matrix)) % logical std
		ix = 1:length(linear_vals); %use them all, pass up the stack
		FP = cumsum(~linear_vals);
		TP = cumsum(linear_vals);
	else % +1 -1 std
		ix = find(linear_vals ~= 0); 
		TP = cumsum(linear_vals(ix) == 1); 
		FP = cumsum(linear_vals(ix) ==-1);
	end
	PR = [TP ./ (TP + FP), TP];
end
function [back_rate, back_num] = compute_background_rate(standard, sga)
	% computes the rate of TP / (TP+FP) for SCREENED pairs
	% Assumes all queries are screened by all arrays
	% Also assumes standard is symmetric

	Q = ismember(standard.orfs, strip_annotation(sga.Cannon.Orf(sga.Cannon.isQuery)));
	A = ismember(standard.orfs, strip_annotation(sga.Cannon.Orf(sga.Cannon.isArray)));
	TP = sum(sum(standard.matrix(Q,A)== 1));
	if islogical(standard.matrix)
		TN = sum(sum(~standard.matrix(Q,A)));
	else
		TN = sum(sum(standard.matrix(Q,A)==-1));
	end

	back_num = TP;
	back_rate = TP / (TP+TN);

	% spit out recall summaries


end
