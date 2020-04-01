function[map] = Hash(map, cell_array)
%function[map] = Hash(map[], cell_array)
% if cell_array has one column, contents are keys and indecies are values
% if cell_array has two columns, col1 are keys and col2 are values
% more than 2 columns are invalid
% (if rows < cols, Hash will transpose)

if(size(cell_array,1) < size(cell_array,2))
	cell_array = cell_array';
end

assert(size(cell_array,2) <= 2);

if(isempty(map))
	map = java.util.HashMap(length(cell_array));
end

if(size(cell_array,2) == 1)
	for i=1:length(cell_array)
		map.put(java.lang.String(cell_array{i}), java.lang.Integer(i));
	end
elseif(size(cell_array,2) == 2)
	for i=1:length(cell_array)
		map.put(java.lang.String(cell_array{i,1}), java.lang.String(cell_array{i,2}));
	end
end
