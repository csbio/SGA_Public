function[list, suffix] = StripOrfs(list, N)
%[list, suffix] = StripOrfs(list, ['first' || 'last'])
% removes _xxx from orf names
% returns the suffix if you ask
% first (default): AA_bb_cc -> AA
% last           : AA_bb_cc -> AA_bb


   if ~exist('N', 'var')
       N = 'first';
	end

	if ~iscell(list)
		[list, suffix] = strip_annotation(list, N);
	else
		suffix = cell(size(list));
		for i=1:numel(list)
			[list{i}, suffix{i}] = strip_annotation(list{i}, N);
		end
	end
end

function[short_orf, suffix] = strip_annotation(orf_string, N)
% from SGA/Main/extra
% Strips off the LAST annotation by default (eg. Orf_sn_rep -> Orf_sn)

    ix = strfind(orf_string, '_');
    if isempty(ix)
        short_orf = orf_string;
        suffix = '';
        return
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

