function linked_orfs = get_linked_ORFs(orfname,genenames,bp_dist,chrom_coords,lfid)

    % Print the name and path of this script
    p = mfilename('fullpath');
    log_printf(lfid, '\nGet linked ORFs script: %s\n\n',p);
    
    linked_orfs = [];

    currgene = strmatch(orfname,genenames,'exact');
    if isempty(currgene)
        return;
    end

    if length(currgene) > 1
        pause; 
    end
    
    % Find all the genes on this chromosome
    indt = find(chrom_coords(:,1) == chrom_coords(currgene,1));
    
    % Sort coordinates by chromosome & position
    [t, chrom_order] = sortrows(chrom_coords, [1 2 3]);

    % Find linked ORFs
    tkb = find(abs(chrom_coords(indt,2)-chrom_coords(currgene,2)) < bp_dist);

    if ~isempty(tkb)
        
        [v,ind] = min(chrom_coords(indt(tkb),2)-chrom_coords(currgene,2));
        lstart = find(chrom_order == indt(tkb(ind)));
        
        [v,ind] = max(chrom_coords(indt(tkb),2)-chrom_coords(currgene,2));
        lend = find(chrom_order == indt(tkb(ind)));
        
        linked_orfs = genenames(chrom_order(lstart:lend));
        
    end



