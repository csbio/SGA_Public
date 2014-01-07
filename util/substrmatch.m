function[indicies, bools] = substrmatch(str, cellary)
%function[vec] = substrmatch(str, cellary)
% see also: Code/cellgrep

bools = ~cellfun(@isempty, strfind(cellary, str));
indicies = find(bools);

end
