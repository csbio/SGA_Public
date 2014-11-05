function[map] = Hash(map, cell_array)
%function[map] = Hash(map[], cell_array)
% when adding to an existing map, indecies start over at i=1

if(isempty(map))
	map = java.util.HashMap(length(cell_array));
end

for i=1:length(cell_array)
	map.put(java.lang.String(cell_array{i}), java.lang.Integer(i));
end
