function[HashMap] = hash_strings(cellarr)
%function[HashMap] = hash_strings(cellarr)
% returns a hash map with strings hashed to their index in cellarr

HashMap = java.util.HashMap(length(cellarr));
for i=1:length(cellarr)
	HashMap.put(java.lang.String(cellarr{i}), java.lang.Integer(i));
end
