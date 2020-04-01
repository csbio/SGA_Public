function [ids] = apply_map(map, strings, table)
%function [ids] = apply_map(map, strings, table)
% returns IDS as defined by map,
% if "table" is provided, returns values 
% with NaNs inserted for id=0 (not in map...)
	ids = zeros(size(strings));
	for i=1:numel(strings)
		ix = map.get(strings{i});
		if ~isempty(ix)
			ids(i) = ix;
		end
	end

	% replace ID return value with lookup values
	% assuming 
	if exist('table', 'var')
		assert(iscell(table))
		vals = ids>0;
		ret = cell(size(ids));
		ret(vals) = table(ids(vals));
		% fill in 
		ret(~vals) = {''};
		ids = ret;
	end
end