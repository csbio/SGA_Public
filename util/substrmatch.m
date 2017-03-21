function[indicies, bools] = substrmatch(str, cellary)
% function[indicies, bools] = substrmatch(str, cellary)
% see also: Code/cellgrep

bools = ~cellfun(@isempty, strfind(cellary, str));
indicies = find(bools);

end
