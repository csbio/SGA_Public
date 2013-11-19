function string = join_by_delimiter(pieces, delimiter)
%function string = join_by_delimiter(pieces, delimiter)

if isempty(pieces)  % no pieces to join, return empty string
    string = '';
else  
    string = pieces{1};
end

l = length(pieces);
p = 1;

while p < l
    p = p+1;
    string = [string delimiter pieces{p}];
end