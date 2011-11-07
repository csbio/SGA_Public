function inds = multistrmatch(s1, s2,exactness,keepone,verbose)

inds = [];
for i = 1 : length(s1)
    if exactness == 1
        x = strmatch(s1{i}, s2,'exact');
    else
        x = strmatch(s1{i}, s2);
    end
    if isempty(x)
        inds = [inds;0];
        if verbose, fprintf('Could not find %s, index set to 0\n', s1{i}); end
    elseif length(x) == 1
        inds = [inds;x];
    else
        if keepone == 0
            if verbose, fprintf('Multiple occurences for %s, saved all\n', s1{i}); end
            inds = [inds; x];
        else
            inds = [inds; x(1)];
            if verbose, fprintf('Multiple occurences for %s, saved the first one\n', s1{i}); end
        end
    end
end