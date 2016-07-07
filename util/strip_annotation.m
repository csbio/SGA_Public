function[heads, tails] = split_annotation(strain_list, pos)
% function[heads, tails] = split_annotation(strain_list, [pos='first'])
%
% Splits strings in the strain list at '_' and returns the pieces.
% strain_list may also be a single string to operate on.
%
% In the case of multiple occurrences of '_' you can specify the break position with asecond arg.
% pos may be a string in {'first', 'last'} default: 'first'
% pos may also be an integer indicating which '_' to break at
% if pos > the number of '_', it will fall back to 'last'

   if ~exist('pos', 'var')
       pos = 'first';
	end

	if ~iscell(strain_list)
		[heads, tails] = single_string(strain_list, pos);
	else
		heads = cell(size(strain_list));
		tails = cell(size(strain_list));
		for i=1:numel(strain_list)
			[heads{i}, tails{i}] = single_string(strain_list{i}, pos);
		end
	end
end

function[short_orf, suffix] = single_string(orf_string, N)

    ix = strfind(orf_string, '_');
    if isempty(ix)
        short_orf = orf_string;
        suffix = '';
    else
        if strcmp('first', N)
            N = 1;

        elseif strcmp('last', N)
            N = length(ix);

        elseif N>length(ix)
            N = length(ix);
        end
        short_orf = orf_string(1:ix(N)-1);
        suffix = orf_string(ix(N)+1:end);
    end
end

