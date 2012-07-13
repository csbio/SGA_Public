function[map] = Hash(map, cell_array)
%function[map] = Hash(map[], cell_array)

if(isempty(map))
	map = java.util.HashMap(length(cell_array));
end

for i=1:length(cell_array)
	map.put(java.lang.String(cell_array{i}), java.lang.Integer(i));
end
