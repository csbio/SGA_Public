function[short_orf] = strip_annotation(orf_string)
%function[short_orf] = strip_annotation(orf_string)
% splits off the LAST annotation (eg. Orf_sn_rep -> Orf_sn)
    ix = strfind(orf_string, '_');
    if(isempty(ix))
        short_orf = orf_string;
    else
        short_orf = orf_string(1:ix(end)-1);
    end
end
