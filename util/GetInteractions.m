function [orfs, coms, eps, pvl] = GetInteractions(sga, query, int_type)
%function [orfs, coms, eps, pvl] = GetInteractions(sga, query, int_type)
% assumes the query exists, and is named exactly
% TODO respect isArray

	if ~exist('int_type', 'var')
		int_type = 'neg';
	end

	THRESH = 0.08;

	assert(ismember(int_type, {'neg'})); % TODO pos, both

	ixq = sga.Cannon.Map.get(query);
	if(isempty(ixq) || ~strcmp(sga.Cannon.Orf{ixq}, query))
		fprintf('warning Cannon.Map inconsistency... doing slow lookup\n')
		ixq = unique([strmatch(query, sga.Cannon.Common); ...
			strmatch(query, sga.Cannon.Common)])
	end

	if strcmp(int_type, 'neg')
		ixa = sga.eps(ixq,:) < -THRESH & sga.pvl(ixq,:) < 0.05;
	else
		fprintf('not implemented yet\n');
		return
	end

	orfs = sga.Cannon.Orf(ixa);
	coms = sga.Cannon.Common(ixa);
	eps = sga.eps(ixq,ixa);
	pvl = sga.pvl(ixq,ixa);
end
