function[sga] = score_trigenic_interactions(sga, assignments)
%function[sga] = score_trigenic_interactions(sga, pairing_file)
% Replaces the scores in a standard sga struct with quant trigenics
% and unsets queries acording to filter rules.
%  removed automatching scheme (now uses find_players)
%  and updated to include data from more inclusive fitness standard
%    scoring model is based on 
%    ExtractQuantitativeTrigenicInteractions_10[_withP].m
%    Model E_abc = E_xx - F_a*E_bc - F_b*E_ac

	fitness = [strip_annotation(sga.Cannon.Orf) num2cell(sga.fit)];

	% Step one, annotate the struct, we don't want to do this twice to the same data
	if(isfield(sga, 'TriGenScored'))
		fprintf('This structure has been scored for trigenics. Skipping...\n');
		return
	else
		sga.TriGenScored = date();
	end

	% We need a map without strain IDS
	NameMap = Hash([], strip_annotation(sga.Cannon.Orf));
	FitNameMap = Hash([], fitness(:,1));

	% Gather the DM_queries and Controls
   % each of the first three return arguments are Nx3
   % orfs and common names are strings (cell)
   % and ix_all are absolute indexes into sga.Cannon, sga.eps, sga.pvl etc
   % first column is double mutant, then singleA then singleB
	[orf_all, com_all, ix_all, sga] = find_trigenic_players(sga, assignments);

	% If these are not unique we may have hash mapping problems, check BEFORE SET OPS
	if length(unique(orf_all(:,1))) < length(orf_all(:,1))
		fprintf('warning, there are some non-unique DM queries\n')
		keyboard
	end

	% find_players knows to toss incomplete sets, we must here unset any dm_query that didn't come back
	tms_ix = substrmatch('+', sga.Cannon.Orf); % all dms
	tms_ix(substrmatch('YDL227C', sga.Cannon.Orf(tms_ix))) = []; % remove controls
	tms_ix(ismember(sga.Cannon.Orf(tms_ix), orf_all(:,1))) = []; % remove found dm_qs
	if(length(tms_ix) > 0)
		fprintf('masking %d dm_queries for incomplete sets\n', length(tms_ix));
		sga.Cannon.isQuery(tms_ix) = false; % unset the remaining dm_qs
		keyboard
	end

	% score complete sets
   % trgenic scoring equation:
   % double mutant query profile (vector of triple "epsilons") minus:
   %     single mutant A profile, scaled by smf of B   minus
   %     single mutant B profile, scaled by smf of A
   for i=1:size(ix_all,1)
      [A_fit, B_fit] = GetSingleFitness(orf_all(i,:), fitness, NameMap);
      sga.eps(ix_all(i,1),:) = sga.eps(ix_all(i,1),:) - B_fit*sga.eps(ix_all(i,2),:) - A_fit*sga.eps(ix_all(i,3),:);
   end
end

function[A_fit, B_fit] = GetSingleFitness(queries, fitness, FitNameMap)
	A_id = FitNameMap.get(queries{2});
	B_id = FitNameMap.get(queries{3});

	% logic to map sm_control strains that weren't specifically assayed for fitness
	% to corresponding "normal" sm strains is included in the script that builds
	% the fitness standard. Therefore, if the fitness is not already in the file,
	% it doesn't exist. 

	if(isempty(A_id))
		A_fit = 1;
	else
		A_fit = fitness{A_id,2};
	end

	if(isempty(B_id))
		B_fit = 1;
	else
		B_fit = fitness{B_id,2};
	end
end
