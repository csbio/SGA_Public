function tokens = split_by_delimiter(string, delimiter)
% function tokens = split_by_delimiter(cell_str, delimiter)
% splits a string or a cell array of strings by delimiter
% 'tokens' will be a cell array with 
% rows=length(string)
% cols= max(# of delims) +1

   if ~iscell(string) 
      string = {string};
   end

   tokens = {};
   for i = 1:length(string)
      remain = string{i};
      tokens{i,1} = '';
      
      j = 1;
      while true
         [str, remain] = strtok(remain, delimiter);
         if ~isempty(str)
            tokens{i,j} = str;
            j = j+1;
         else
            break;
         end
      end
   end
end

