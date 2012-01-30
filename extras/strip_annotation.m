function[short_orf] = strip_annotation(orf_string, N)
%function[short_orf] = strip_annotation(orf_string, [N={'first', 'last', #N})
%
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
