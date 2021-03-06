function[result] = OrfToCommon(cellarr, mappingfile)
%function[result] = OrfToCommon(cellarr, [mappingfile | 'quiet'])
% if the inputs have _tails, preserve for the output

% parse the mapping file
if(~exist('mappingfile', 'var'))
	% use defalut, with a warning
	mappingfile = [get_SGAROOT() '/refdata/strain_orf_common_allele_160425.txt'];
	fprintf('using default mapping file: %s\n', mappingfile);
elseif(strcmp(mappingfile, 'quiet'))
	% use defalut, but no warning
	mappingfile = [get_SGAROOT() '/refdata/strain_orf_common_allele_160425.txt'];
end



fid = fopen(mappingfile, 'r');
A = textscan(fid, '%s%s%s%s', 'Delimiter', '\t', 'ReturnOnError', false);
fclose(fid);

% handle tails
[cellarr, tails] = strip_annotation(cellarr, 'first');


result = cell(size(cellarr));
cellarr = cellarr(:);
map = hash_strings(A{2});

for i=1:length(cellarr)
	ix = map.get(cellarr{i});
	if(isempty(ix) || isempty(A{3}{ix}))
		% Check for possible double query
		ix_plus = strfind(cellarr{i}, '+');
		if isempty(ix_plus)
			result{i} = cellarr{i};
		elseif(length(ix_plus) == 1)
			ix1 =map.get(cellarr{i}(1:ix_plus-1));
			if(isempty(ix1) || isempty(A{3}{ix1}))
				r1 = cellarr{i}(1:ix_plus-1);
			else
				r1 = A{3}{ix1};
			end

			ix2 =map.get(cellarr{i}(ix_plus+1:end));
			if(isempty(ix2) || isempty(A{3}{ix2}))
				r2 = cellarr{i}(ix_plus+1:end);
			else
				r2 = A{3}{ix2};
			end
			result{i} = [r1 '+' r2];
		else
			error('warning, too many +\n');
		end
	else
		result{i} = A{3}{ix};
	end
end


% replace tails
for i=1:length(result)
	if ~isempty(tails{i})
		result{i} = [result{i} '_' tails{i}];
	end
end
