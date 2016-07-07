function[] = PRCurves_from_files(file_list, input_label, standard, std_label)
%function[] = PRCurves_from_files(file_list, input_label, standard, std_label)

	lists = cell(size(file_list));
	values= cell(size(file_list));
	for i=1:length(lists)
		fid = fopen(file_list{i}, 'r');
		% Determine number of columns (2 and 3 supported, tabs)
		A = textscan(fgetl(fid), '%s', 'Delimiter', '\t', 'ReturnOnError', false);
		fseek(fid,0,-1); % rewind line	
		COLS = length(A{1});
		if(COLS == 2) 
			A = textscan(fid, '%s%s', 'Delimiter', '\t', 'ReturnOnError', false, 'BufSize', 2^20); 
		elseif(COLS == 3);% We've got a numerical column, sort it
			A = textscan(fid, '%s%s%f', 'Delimiter', '\t', 'ReturnOnError', false, 'BufSize', 2^20); 
			if ~issorted(A{3})
				fprintf('sorting file %d\n', i);
				[A{3}, ix] = sort(A{3}, 'descend');
				A{1} = A{1}(ix);
				A{2} = A{2}(ix);
			end
			values{i} = A{3};
		else
			error('unsupported number of columns or wrong delimiter. Make me smarter.\n');
		end
		fclose(fid);

		% now we only need cols 1 and 2
		lists{i} = [A{1} A{2}];
	end


	% (Optional) Restrict PR Evaluation to the top (min(length))
	% better to do this before calling this script but this is better than nothing for head-to-head
	%{
	num_pairs = min(cellfun(@length, lists));
	for i=1:length(lists)
		lists{i} = lists{i}(1:num_pairs,:);
	end
	%}

	% (Optional) Convert Alleles to Orfs
	% This code originally written for Allele-Allele similarity standard evaluation.
	for i=1:length(lists)
		%lists{i} = [AlleleToOrf(lists{i}(:,1)) AlleleToOrf(lists{i}(:,2))];
		lists{i} = [StripOrfs(AlleleToOrf(lists{i}(:,1))) ...
		            StripOrfs(AlleleToOrf(lists{i}(:,2)))];
	end

	% Remove pairs that contain an ORF not in our std.
	% they might not be the same length after this.
	for i=1:length(lists)
		keep = ismember(lists{i}(:,1), standard.orfs) & ismember(lists{i}(:,2), standard.orfs);
		lists{i} = lists{i}(keep,:);
	end

	% Evaluate the curves
	PRs = cell(size(lists));
	for i=1:length(lists)
		PRs{i} = help_eval_ppi_pairs(standard, lists{i});
	end

	% Draw the curves
	CM = colormap();
	for i=1:length(lists)
		colors{i} = CM(floor(size(CM,1)/(length(lists)+1))*i,:);
		semilogx(PRs{i}(:,2), PRs{i}(:,1), 'Color', colors{i}, 'LineWidth', 2);
		hold on;
	end

	xlabel('Recall');
	ylabel('Precision');
	title(sprintf('%s predicts %s', input_label, std_label), 'Interpreter', 'None');
	legend(file_list, 'Interpreter', 'None', 'Location', 'SW');
	hold off

end



% sub-function lifted from compare_sga_structs on Apr 11, 2013
function[PR, ix] = help_eval_ppi_pairs(standard, pair_list)

   name_map = hash_strings(standard.orfs);
   
   rows = zeros(size(pair_list, 1),1);
   cols = zeros(size(pair_list, 1),1);
   % none of these should return empty...
   for i=1:size(pair_list, 1); 
      rows(i) = name_map.get(pair_list{i,1});
      cols(i) = name_map.get(pair_list{i,2});
   end 
   linear_ind = sub2ind(size(standard.matrix), rows, cols);
   linear_vals = standard.matrix(linear_ind);
   ix = find(linear_vals ~= 0); 
   TP = cumsum(linear_vals(ix) == 1); 
   FP = cumsum(linear_vals(ix) ==-1);
   PR = [TP ./ (TP + FP), TP];
end

