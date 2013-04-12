function[genome_ix] = ChromSort(orf_array)
%[ix] = ChromSort(orf_array)
% Sorts a cell array of ORFs into chromosome (nearly alpha) order 
% < L CHR1 R > < L CHR2 R > ...


% Automatically detects DM controls (e.g. ignores YDL227C)
% the return in indicies, so we can mangle the inputs if we need
% true DM queries will get sorted according to position 1
DMC = substrmatch('YDL227C+', orf_array);
[trsh, orf_array(DMC)] = SplitOrfs(orf_array(DMC));

% these will sort fine already
% DMC = substrmatch('+YDL227C', orf_array);
% [orf_array(DMC), trsh] = SplitOrfs(orf_array(DMC));



CHR = 'ABCDEFGHIJKLMNOP';
ptr = 1;

genome_ix = nan(size(orf_array));
for i=1:length(CHR)
	left_ix = strmatch(['Y' CHR(i) 'L'], orf_array);% extract absolute ix
	[~, ix_ix] = sortrows(orf_array(left_ix));	
	left_ix = left_ix(ix_ix(end:-1:1));							   % rearrange absolute ix
	genome_ix(ptr:ptr+length(left_ix)-1) = left_ix;
	ptr = ptr + length(left_ix);


	right_ix = strmatch(['Y' CHR(i) 'R'], orf_array);% extract absolute ix
	[~, ix_ix] = sortrows(orf_array(right_ix));	
	right_ix = right_ix(ix_ix);							   % rearrange absolute ix
	genome_ix(ptr:ptr+length(right_ix)-1) = right_ix;
	ptr = ptr + length(right_ix);
end


other_ix = setdiff(1:length(orf_array), genome_ix);
[~, ix_ix] = sortrows(orf_array(other_ix));
other_ix= other_ix(ix_ix);
genome_ix(ptr:ptr+length(other_ix)-1) = other_ix;
ptr = ptr + length(other_ix);

assert(ptr == length(orf_array)+1);
