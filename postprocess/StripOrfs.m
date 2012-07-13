function[list] = StripOrfs(list, N)
%[list] = StripOrfs(list, ['first' || 'last'])
% removes _xxx from orf names
% first			 : AA_bb_cc -> AA
% last (default): AA_bb_cc -> AA_bb


   if ~exist('N', 'var')
       N = 'last';
	end

	if ~iscell(list)
		list = strip_annotation(list, N);
	else
		for i=1:length(list)
			list{i} = strip_annotation(list{i}, N);
		end
	end

end

function[short_orf] = strip_annotation(orf_string, N)
% from SGA/Main/extra
% Strips off the LAST annotation by default (eg. Orf_sn_rep -> Orf_sn)
% you may also provide explicit instructions 'first', 'last', or N
%
% strip_annotation() will remove the last annotation if given too small or too large an N

    ix = strfind(orf_string, '_');
    if isempty(ix)
        short_orf = orf_string;
        return
    else

        if exist('N', 'var')
            if strcmp('first', N)
                N = 1;

            elseif strcmp('last', N)
                N = length(ix);

            elseif N>length(ix)
                N = length(ix);
            end
        else
            N = length(ix);
        end
        short_orf = orf_string(1:ix(N)-1);
    end
end

