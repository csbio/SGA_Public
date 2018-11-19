function[orf_all, com_all, ix_all, sga] = find_trigenic_players(sga, assignment_file)
%function[orf_all, com_all, ix_all, sga] = find_trigenic_players(sga, assignment_file)
% indicies are 'absolute'
% results are Nx3    DM, paralog1, paralog2
% Will update if isQuery has been restricted, but not expanded.
% To reflect expanded isQuery delete the .assignments field.
% if the assignments struct is present, it gets passed back out.
% tosses any query with insufficient controls
% 
% assignment_file is a file Nx4, of strain ids dm, c1, c2, set/class
% returns SGA with assignment struct params RESET
	% you can skip / force with (..., [], true)

	% look for the assignment file path in the struct if not passed in
	if ~exist('assignment_file', 'var') && isfield(sga, 'assignment_file')
		assignment_file = sga.assignment_file;
	end

	% if we've saved assignments already, just pass them back out
	if isfield(sga, 'assignments')
		orf_all = sga.assignments.orf_all;
		com_all = sga.assignments.com_all;
		ix_all  = sga.assignments.ix_all;

		% isQuery may have changed since assignment mapping, so strip incomplete sets
		% e.g. project set selection 
		valid_queries = ismember(orf_all, sga.Cannon.Orf(sga.Cannon.isQuery));
		complete_sets = all(valid_queries,2);

		orf_all = orf_all(complete_sets,:);
		com_all = com_all(complete_sets,:);
		ix_all = ix_all(complete_sets,:);
		return
	end
	% set fresh from file
	if exist('assignment_file', 'var')
		[~, tails] = strip_annotation(sga.Cannon.Orf);
		% unset the non-query tails to respect isQuery
		tails(~sga.Cannon.isQuery) = {'NA'};
		tailmap = Hash([], tails);

		% load up the mappings
		% dm sm1 sm2 set/class
		fid = fopen(assignment_file, 'r');
		assigns = textscan(fid, '%s%s%s%s', 'Delimiter', '\t', ...
				'HeaderLines', 1, 'ReturnOnError', false);
		assigns = [assigns{:}];
		
		ix_all = apply_map(tailmap, assigns(:,1:3));
		toss = sum(ix_all == 0, 2)>0;
		fprintf('%d sets from assigments not found in data (may include unset queries)\n', sum(toss));
		ix_all = ix_all(~toss, :);
		class_all = assigns(~toss,4);

		orf_all = cell(size(ix_all));
		orf_all(:) = sga.Cannon.Orf(ix_all(:));
		com_all = OrfToCommon(orf_all);

		% make sure each dm strain is listed once, so we don't subtract
		% twice when scoring
		[dup_list, dup_count, dup_ix] = FindDuplicates(orf_all(:,1));
		if ~isempty(dup_list)
			fprintf('promiscuous dm strains!\n');
			fprintf('%d total strains appear %d total times (%d total sets)\n', length(dup_list), sum(dup_count), size(orf_all,1));
			fprintf('"dbcont" to continue, tossing strains with ambiguous sets.\n');
			toss_ix = [dup_ix{:}];
			keyboard

			orf_all(toss_ix,:) = [];
			com_all(toss_ix,:) = [];
			ix_all(toss_ix,:) = [];
			class_all(toss_ix,:) = [];
		end

		sga.assignments = struct();
		sga.assignments.orf_all = orf_all;
		sga.assignments.com_all = com_all;
		sga.assignments.ix_all = ix_all;
		sga.assignments.class_all = class_all;
		sga.assignment_file = assignment_file;

	% no assignment information, try to auto-detect (older style)
	else

		allow_incomplete = true;
		if allow_incomplete
			fprintf('warning find_trigenic_players returning incomplete sets\n');
		end

		% update the Map, just in case
		sga.Cannon.Map = Hash([], sga.Cannon.Orf);

		Qorfs = sga.Cannon.Orf(sga.Cannon.isQuery);
		DM_orfs = Qorfs(substrmatch('+', sga.Cannon.Orf(sga.Cannon.isQuery)));
		DM_case = setdiff(DM_orfs, DM_orfs(substrmatch('YDL227C', DM_orfs)));
		DM_ctrl = setdiff(DM_orfs, DM_case);

		ix_all = zeros(length(DM_case),3); % DM A B
		for i=1:size(ix_all,1)
			ix_all(i,1) = sga.Cannon.Map.get(DM_case{i});
			[A, B] = SplitOrfs(strip_annotation(DM_case{i}), '+');
			A_ctrl = substrmatch([A '+YDL227C'], DM_ctrl);
			B_ctrl = substrmatch(['YDL227C+' B], DM_ctrl);
			if(~isempty(A_ctrl))
				ix_all(i,2) = sga.Cannon.Map.get(DM_ctrl{A_ctrl});
			end
			if(~isempty(B_ctrl))
				ix_all(i,3) = sga.Cannon.Map.get(DM_ctrl{B_ctrl});
			end
		end

		if ~allow_incomplete % overwrite without checking
			% strip out any set that wasn't complete
			toss = sum(ix_all ==0, 2) >0;
			ix_all = ix_all(~toss,:);
			fprintf('tossing %d sets\n', sum(toss));
			fprintf('controls left over: %d\n', length(setdiff(DM_ctrl, sga.Cannon.Orf(ix_all(:)))));
			orf_all = sga.Cannon.Orf(ix_all);
			com_all = sga.Cannon.Common(ix_all);
		else
			% ix_all may have some zeros
			orf_all = cell(size(ix_all));
			com_all = cell(size(ix_all));
			orf_all(ix_all ~= 0) = sga.Cannon.Orf(ix_all(ix_all ~= 0));
			com_all(ix_all ~= 0) = sga.Cannon.Common(ix_all(ix_all ~= 0));
		end


	% above code is keyed on dm_queries, so we still need to find a home for left over controls
		if allow_incomplete 
			% we have no basis for matching controls without triple mutants so these will have no partners
			DM_leftover = setdiff(DM_ctrl, orf_all(~cellfun(@isempty, orf_all)));
			DM_leftover_A = DM_leftover(substrmatch('+YDL227C', DM_leftover));
			DM_leftover_B = DM_leftover(substrmatch('YDL227C+', DM_leftover));

			orf_all_append = cell(length(DM_leftover_A),3);
			com_all_append = cell(length(DM_leftover_A),3);
			ix_all_append = zeros(length(DM_leftover_A),3);
			for i=1:length(DM_leftover_A)
				ix_all_append(i,2) = sga.Cannon.Map.get(DM_leftover_A{i});
				orf_all_append{i,2} = sga.Cannon.Orf{ix_all_append(i,2)};
				com_all_append{i,2} = sga.Cannon.Common{ix_all_append(i,2)};
			end
			ix_all = [ix_all; ix_all_append];
			orf_all = [orf_all; orf_all_append];
			com_all = [com_all; com_all_append];

			ix_all_append = zeros(length(DM_leftover_B),3);
			orf_all_append = cell(length(DM_leftover_B),3);
			com_all_append = cell(length(DM_leftover_B),3);
			for i=1:length(DM_leftover_B)
				ix_all_append(i,3) = sga.Cannon.Map.get(DM_leftover_B{i});
				orf_all_append{i,3} = DM_leftover_B{i};
				com_all_append{i,3} = sga.Cannon.Common{ix_all_append(i,3)};
			end
			ix_all = [ix_all; ix_all_append];
			orf_all = [orf_all; orf_all_append];
			com_all = [com_all; com_all_append];
		end
	end
end
