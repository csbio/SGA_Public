function [] = compare_sga_structs(sga, labels, standard, std_label, varargin)
%function[] = compare_sga_structs(sga, labels, standard, std_label, ['intersect', 'random_allele'])

	PLOTS = true;
	SINGLE_PLOT = false;
	%Colors = {'r', 'b', 'g', 'm', 'k', 'c', 'y'};

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

	% OPTIONAL intersect queries and arrays
	if(ismember('intersect', varargin))
		all_orfs = sga{1}.Cannon.Orf;
		for i=2:length(sga)
			all_orfs = intersect(all_orfs, sga{i}.Cannon.Orf);
		end
		all_orfs = unique(all_orfs);
		for i=1:length(sga)
			sga{i}.Cannon.isQuery(~ismember(sga{i}.Cannon.Orf, all_orfs)) = false;
			sga{i}.Cannon.isArray(~ismember(sga{i}.Cannon.Orf, all_orfs)) = false;
		end
	fprintf('Intersection:\n\tQueries: %d\n\tArrays: %d\n', ...
                 sum(sga{1}.Cannon.isQuery), sum(sga{1}.Cannon.isArray));
	end

	
	for i=1:length(sga)
		% mask out multiple alleles (Q & A) if so asked
		if(ismember('random_allele', varargin))
			sga{i} = reduce_alleles(sga{i});
		end
		%remove nans (not found in std)
		% mask out any queries or arrays not in the standard
		sga{i}.Cannon.isQuery(~ismember(StripOrfs(sga{i}.Cannon.Orf), standard.orfs)) = false;
		sga{i}.Cannon.isArray(~ismember(StripOrfs(sga{i}.Cannon.Orf), standard.orfs)) = false;

		% perform the analysis
      PR_neg = help_curve_the_pr(sga{i}, standard, true);
      PR_pos = help_curve_the_pr(sga{i}, standard, false);

		% draw the results

		if(PLOTS)
			% background, draw this for each SGA in case they are different.
			% background is the same for pos and negative.
			colors{i} = CM(floor(size(CM,1)/(length(sga)+1))*i,:);
			backgrounds(i) = compute_background_rate(standard, sga{i});

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
		% Plot the backgrounds afterward so they don't show up in the legend
		for i=1:length(sga)
			xx = get(gca, 'XLim');
			plot(xx, [backgrounds(i) backgrounds(i)], ...
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
			plot(xx, [backgrounds(i) backgrounds(i)], ...
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
   Qname = StripOrfs(sga.Cannon.Orf(sga.Cannon.isQuery), 'first');
   Aname = StripOrfs(sga.Cannon.Orf(sga.Cannon.isArray), 'first');
   
   [tmp, ix] = sort(abs(v), 'descend');
   r = r(ix);
   c = c(ix);
   v = v(ix);

   [PR, ix] = help_eval_ppi_pairs(standard, [Qname(r) Aname(c)]);
   PR = [v(ix) PR];

end
function [PR, ix] = help_eval_ppi_pairs(standard, pair_list)

	name_map = Hash(java.util.HashMap(), standard.orfs);
   
	rows = zeros(size(pair_list, 1),1);
	cols = zeros(size(pair_list, 1),1);
	% none of these should return empty...
	for i=1:size(pair_list, 1); 
		rows(i) = name_map.get(pair_list{i,1});
		cols(i) = name_map.get(pair_list{i,2});
	end 
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
function [back] = compute_background_rate(standard, sga)
	% computes the rate of TP / (TP+FP) for SCREENED pairs
	% Assumes all queries are screened by all arrays
	% Also assumes standard is symmetric

	Q = ismember(standard.orfs, StripOrfs(sga.Cannon.Orf(sga.Cannon.isQuery)));
	A = ismember(standard.orfs, StripOrfs(sga.Cannon.Orf(sga.Cannon.isArray)));
	TP = sum(sum(standard.matrix(Q,A)== 1));
	if islogical(standard.matrix)
		TN = sum(sum(~standard.matrix(Q,A)));
	else
		TN = sum(sum(standard.matrix(Q,A)==-1));
	end
	back = TP / (TP+TN);

	% ? what p
end
