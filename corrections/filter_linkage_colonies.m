%%
% FILTER_LINKAGE_COLONIES - generates the list of colonies affected by linkage
%
% Inputs:
%   sgadata - structure containing all the colony size data
%   linkagefile - name of the file containing the coordinates of the linkage windows
%               <ORF | Strain> <CHRM> <START> <END>
%   coord_file - name of the file containing coordinates to map arrays
%               <ORF> <CHRM> <START> <END>
%   all_arrays, all_querys - precompiled integer IDs from compute_sgacore.m
%   query_map, array_map - precomputed maps of indecies into sgadata
%   wild_type - string ID of the wild-type strain
%
% Outputs:
%   all_linkage_cols - indices of the colonies affected by linkage
%
% Authors: Chad Myers (cmyers@cs.umn.edu), 
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca),
%          Benjamin VanderSluis (bvander@cs.umn.edu)
%
% Linkage matching behavior
%    First find global linkage colonies (e.g. linked to query markers)
%    Then for each query:
%        a) Look for strain match in linkage file (may be > 1 entry)
%        b) If that's not found, look for defined windows for each ORF 
%        c) If that's not found, map a default window for each ORF
%
%        For b and c iterate over N orfs, delimited by '+' in the strainid
%%

function [all_linkage_cols, non_spec]  = filter_linkage_colonies_multiregion(sgadata, linkagefile, coord_file, all_querys, all_arrays, query_map, array_map, wild_type, lfid)

    % Print the name and path of this script
    p = mfilename('fullpath');
    log_printf(lfid, '\nLinkage filter script:\n\t%s\n\n',p);

    match_code_counts = zeros(1,4);
    match_code_labels = {'query strain based linkages:', 'query orf based linkages:', ...
                        'query window based linkages:', 'query linkage failures:'};

    wild_type_id = find(strncmp(wild_type, sgadata.orfnames, length(wild_type))); % partial match 
    all_linkage_cols_bool = false(size(sgadata.arrays)); 

    % ---------------------------------
    % Load chromosomal coordinates 
    % ---------------------------------
    % Format <string ORF> <int CHROMOSOME> <START position> <END position>
    % These are canonical coordinates from SGD, therefore START is not strictly < END
    coord_fid = fopen(coord_file, 'r');
    coord_data = textscan(coord_fid, '%s%d%d%d', 'Delimiter', '\t', 'ReturnOnError', false);
    fclose(coord_fid);
    SGD_coord = struct();
    SGD_coord.orfs = coord_data{1};
    SGD_coord.locations = [coord_data{3} coord_data{4}];
    SGD_coord.locations = sort(SGD_coord.locations, 2, 'ascend'); % enforce low to high assumption
    SGD_coord.locations = [coord_data{2} SGD_coord.locations];     % prepend chrom number

    % ---------------------------------
    % load predefined linkage
    % ---------------------------------
    % Format <string ORF | STRAIN> <int CHROMOSOME> <START position> <END position>
    % START should be LESS THAN END, but this is enforced, not assumed
    linkage_fid = fopen(linkagefile, 'r');
    predef_lnkg_data = textscan(linkage_fid, '%s%d%d%d', 'Delimiter', '\t', 'ReturnOnError', false); 
    fclose(linkage_fid);
    % textscan won't complain if we ask for more columns than actually exist, (just returns 0's)
    % so we need to test for the older, (sans-chromosome) linkage file format explicity
    assert(min(predef_lnkg_data{4}) > 0, 'Check linkgefile format (should be 4 col)');
    predef_lnkg = struct();
    predef_lnkg.orf = predef_lnkg_data{1};
    predef_lnkg.chrom = predef_lnkg_data{2};
    predef_lnkg.coord = [predef_lnkg_data{3} predef_lnkg_data{4}];
    predef_lnkg.coord = sort(predef_lnkg.coord, 2, 'ascend'); % enforce low to high assumption

    % ---------------------------------
    % map ARRAY coordinates 
    % ---------------------------------
    array_orfs_found = 0;
    array_orfs_not_found = 0;
    array_coord = nan(length(all_arrays),3); % CHR START END
    for i = 1:length(all_arrays)
        array_orf = strip_annotation(sgadata.orfnames{all_arrays(i)});
        ix = strcmp(array_orf, SGD_coord.orfs);

        if any(ix)
            array_orfs_found = array_orfs_found + 1;
            array_coord(i,:) = SGD_coord.locations(ix,:);
        else
            array_orfs_not_found = array_orfs_not_found + 1;
        end
    end

    log_printf(lfid, 'Array coordinates mapped. \n');
    log_printf(lfid, '\tarray orfs matched:   %d\n', array_orfs_found);
    log_printf(lfid, '\tarray orfs not found: %d\n', array_orfs_not_found);

    % ----------------------------------------------------------
    % global linkage considerations, e.g. marker-linked arrays
    % ----------------------------------------------------------
    
    % Remove lyp and can linkage for all queries
    % (this set of arrays will be empty on the mini1200 array)
    [lyp1_linkage_arrays, ~] = get_linked_arrays('YNL268W', predef_lnkg, SGD_coord, array_coord);
    [can1_linkage_arrays, ~] = get_linked_arrays('YEL063C', predef_lnkg, SGD_coord, array_coord);
    lyp_can_linkage = unique([lyp1_linkage_arrays; can1_linkage_arrays]);
    for i=1:length(lyp_can_linkage)
        all_linkage_cols_bool(array_map{all_arrays(lyp_can_linkage(i))}) = true;
    end 

    % remove ura3 linkage for wild-type
    [ura3_linkage_arrays, ~] = get_linked_arrays('YEL021W', predef_lnkg, SGD_coord, array_coord);
    for i=1:length(ura3_linkage_arrays)
        for j=1:length(wild_type_id)
            all_linkage_cols_bool(intersect(array_map{all_arrays(ura3_linkage_arrays(i))}, ...
                                            query_map{wild_type_id(j)})) = true;
        end
    end

    %if every query has a +, this is triple mutant triple array, and we need to unlink HO globally
    % we also need to unlink URA3 globally, not just for WTS
    if(all(~cellfun(@isempty, strfind(sgadata.orfnames(all_querys), '+'))));
        log_printf(lfid, 'All Queries are DM, removing HO/URA3 globally...\n');
        [ho_linkage_arrays, ~]   = get_linked_arrays('YDL227C', predef_lnkg, SGD_coord, array_coord);
        ura3_ho_linkage = unique([ura3_linkage_arrays; ho_linkage_arrays]);
        for i=1:length(lyp_can_linkage)
            all_linkage_cols_bool(array_map{all_arrays(ura3_ho_linkage(i))}) = true;
        end 
    end

    % pass out all of the global (non-query-specific) linkage colonies
    non_spec = find(all_linkage_cols_bool);

    % -------------------------------------------------------------------------
    % determine specific linkage for each query strain
    % -------------------------------------------------------------------------
    log_printf(lfid, ['Mapping query-specific linkage...\n|', blanks(50), '|\n|']);
    for i = 1:length(all_querys)

        % get all arrays for this strain
        query_string = sgadata.orfnames{all_querys(i)};
        [query_linked_arrays, match_code] = get_linked_arrays(query_string, ...
                                predef_lnkg, SGD_coord, array_coord);

        % convert to global ix and record
        for j = 1:length(query_linked_arrays)
            all_linkage_cols_bool(intersect(array_map{all_arrays(query_linked_arrays(j))}, ...
                                            query_map{all_querys(i)})) = true;
        end

        % Print progress
        match_code_counts(match_code) = match_code_counts(match_code) + 1;
        print_progress(lfid, length(all_querys), i);
    end

    log_printf(lfid, '|\nLinkage Query Process Report\n');
    for i=1:length(match_code_labels)
        log_printf(lfid, '\t%s\t%d\n', match_code_labels{i}, match_code_counts(i));
    end

    all_linkage_cols = find(all_linkage_cols_bool);
end

function [linked_arrays, match_code]  = get_linked_arrays(query_string, predef_lnkg, SGD_coord, array_coord)
% returns indicies into all_arrays
    % query_string will be either an entire strain ID or a single ORF
    % match_code keeps track of how we determined the linkage resion for this strain/orf 
    % 1: we found this entire query string (e.g. strain) in the linkage file
    % 2: we found the stripped version of the string in the linkage file
    % 3: we built a window from coordinates we found for this orf
    % 4: we failed to determine any linkage at all or are still looking
    % 0: internal value, means still looking...
    % double orfs will just get a match code from the second entry

    linked_arrays = []; % WTs at least will match nothing...
    match_code = 4;
    windows = [];

    % First look for exact windows for this strain (skip if this is an orf)
    if any(query_string == '_') % input is a strain
        ix = strcmp(query_string, predef_lnkg.orf);
        if any(ix) % we found this strain
            match_code = 1;
            windows = [predef_lnkg.chrom(ix) predef_lnkg.coord(ix,:)];
        end
    end

    % next we look for the ORF(s) in the predef file
    % even if we have found a strain match, continue looking for an ORF match
    if match_code == 4 or match_code == 1
        orf_list = split_by_delimiter('+', strip_annotation(query_string, 'first'));
        for i=1:length(orf_list)
            this_orf_chrom = upper(orf_list{i}(2)) - '@'; 
            ix = strcmp(orf_list{i}, predef_lnkg.orf);
            if any(ix)
                % we found entries for this orf
                match_code = 2;
                windows = [windows; predef_lnkg.chrom(ix) predef_lnkg.coord(ix,:)];
            else
                % try to map a window to this orf
                [lw, rw, match_code] = lnkg_window_from_coord(orf_list{i}, SGD_coord);
                if match_code == 3
                    windows = [windows; [this_orf_chrom lw rw]];
                end
            end
        end
    end

    % if we have something by now
    if match_code < 4
        linked_arrays = false(size(array_coord,1),1);
        for i=1:size(windows,1)
            chrom = windows(i,1);
            left_win = windows(i,2);
            right_win = windows(i,3);

            linked_arrays((array_coord(:,1) == chrom) & ...
                           ((array_coord(:,2) > left_win & array_coord(:,2) < right_win) | ...
                            (array_coord(:,3) > left_win & array_coord(:,3) < right_win))) = true;
        end
        linked_arrays = find(linked_arrays); % convert to ix
    end
end

function[left_win, right_win, match_code] = lnkg_window_from_coord(orf_string, SGD_coord)
    linkage_dist = 200e3;   % default linkage distance
    left_win = 1; % default values for error
    right_win = 1;
    match_code = 4;
    ix = strcmp(orf_string, SGD_coord.orfs);
    if any(ix)
        match_code = 3;
        left_win = max(SGD_coord.locations(ix,2) - linkage_dist, 1);
        right_win = SGD_coord.locations(ix,3) + linkage_dist;
    end
end
