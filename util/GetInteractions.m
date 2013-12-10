function [orfs, coms, eps, pvl] = GetInteractions(sga, query, int_type, THRESH)
%function [orfs, coms, eps, pvl] = GetInteractions(sga, query, int_type, THRESH)
% assumes the query exists, and is named exactly
% THRESH should be positive (abs())
% type defaults to 'neg' and THRESH to '0.08'

	if ~exist('int_type', 'var')
		int_type = 'neg';
	end

	if ~exist('THRESH', 'var')
		THRESH = 0.08;
	end

	assert(ismember(int_type, {'neg'})); % TODO pos, both

	ixq = sga.Cannon.Map.get(query);
	if(isempty(ixq) || ~strcmp(sga.Cannon.Orf{ixq}, query))
		fprintf('warning Cannon.Map inconsistency... doing slow lookup\n')
		ixq = unique([strmatch(query, sga.Cannon.Common); ...
			strmatch(query, sga.Cannon.Common)])
	end

	if strcmp(int_type, 'neg')
		ixa = sga.eps(ixq,sga.Cannon.isArray) < -THRESH & sga.pvl(ixq,sga.Cannon.isArray) < 0.05;
	else
		fprintf('not implemented yet\n');
		return
	end

	arrays = find(sga.Cannon.isArray)
	orfs = sga.Cannon.Orf(arrays(ixa));
	coms = sga.Cannon.Common(arrays(ixa));
	eps = sga.eps(ixq,arrays(ixa));
	pvl = sga.pvl(ixq,arrays(ixa));
end
