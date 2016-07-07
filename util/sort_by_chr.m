function[sorted_list, sorted_ix] = sort_by_chr(orf_array)
% function[sorted_list, sorted_ix] = sort_by_chr(orf_array)
% Sorts a cell array of ORFs into chromosome (nearly alpha) order 
% < L CHR1 R > < L CHR2 R > ...
%
% Automatically detects DM controls:
% e.g. YDL227C+YPL003W_tm2902 == YPL003W_tm2902 

   % replace these with second orf first, so they sort properly
   orf_array_save = orf_array;
   DMC = substrmatch('YDL227C+', orf_array);
   [~, orf_array(DMC)] = SplitOrfs(orf_array(DMC));

   % preallocate answer, ptr to keep place
   sorted_ix = NaN(size(orf_array));
   ptr = 1;

   CHR = 'ABCDEFGHIJKLMNOP';
   for i=1:length(CHR)
      % left_ix is original indices for left arm of chr i
      left_ix = strmatch(['Y' CHR(i) 'L'], orf_array);    

      % sort them relative to one another
      [~, ix_ix] = sortrows(orf_array(left_ix));	
      left_ix = left_ix(ix_ix(end:-1:1));  % left arm reversed

      % put them into the answer and advance ptr
      sorted_ix(ptr:ptr+length(left_ix)-1) = left_ix;
      ptr = ptr + length(left_ix);


      % repeat for right
      right_ix = strmatch(['Y' CHR(i) 'R'], orf_array);% extract absolute ix
      [~, ix_ix] = sortrows(orf_array(right_ix));	
      right_ix = right_ix(ix_ix); % right arm NOT reversed
      sorted_ix(ptr:ptr+length(right_ix)-1) = right_ix;
      ptr = ptr + length(right_ix);
   end

   % this catches any cases that didn't match or expectations ORF strings
   other_ix = setdiff(1:length(orf_array), sorted_ix);
   [~, ix_ix] = sortrows(orf_array(other_ix));
   other_ix= other_ix(ix_ix);
   sorted_ix(ptr:ptr+length(other_ix)-1) = other_ix;
   ptr = ptr + length(other_ix);

   sorted_list = orf_array_save(sorted_ix);
end
