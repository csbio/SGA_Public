function [map] = hash_strings(cellarr, map)
% function [map] = hash_strings(cellarr, map)
% returns a hash map with strings hashed to their index in cellarr
% if you provide an existing HashMap, keys will be added there instead

   if ~exist('map', 'var')
      map = java.util.HashMap(length(cellarr));
   end
   for i=1:length(cellarr)
      map.put(java.lang.String(cellarr{i}), java.lang.Integer(i));
   end
end
