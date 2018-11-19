function[a, b] = SplitOrfs(cellarr, delim)
%[a, b] = SplitOrfs(cellarr, [delim = '+'])
   % splits on first +, returns pieces

   if ~exist('delim', 'var')
      delim = '+';
   end

   % support a single string also
   if isstr(cellarr)
      [a, b] = single_split(cellarr, delim);
   else
      a = cell(size(cellarr));
      b = cell(size(cellarr));
      for i=1:length(cellarr)
         [a{i}, b{i}] = single_split(cellarr{i}, delim);
      end
   end
end

function[a, b] = single_split(string, delim)
   ix = strfind(string, delim);
   if(isempty(ix))
      a = string;
      b = '';
   else
      a = string(1:ix-1);
      b = string(ix+1:end);
   end
end

