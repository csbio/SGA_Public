function [sorted_list,ix] = sort_by_chr(list)

chr_order = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'};

sorted_list = [];
ix = [];
for c = 1 : length(chr_order)
    indL = strmatch(['Y', chr_order{c}, 'L'], list);
    if ~isempty(indL)
        [listL,ixL] = sortrows(list(indL), -1);
        indL = indL(ixL);
    end
    
    indR = strmatch(['Y', chr_order{c}, 'R'], list);
    if ~isempty(indR)
        [listR,ixR] = sortrows(list(indR));
        indR = indR(ixR);
    end
    
    sorted_list = [sorted_list; list(indL); list(indR)];
    ix = [ix; indL; indR];
      
end
    
