%%
% FILTER_ALL_LINKAGE_COLONIES_QUERYSPECIFIC - generates the list of colonies affected by linkage
%
% Inputs:
%   sgadata - structure containing all the colony size data
%   linkagefile - the name of the file containing the coordinates of the linkage windows
%
% Outputs:
%   all_linkage_cols - indices of the colonies affected by linkage
%
% Authors: Chad Myers (cmyers@cs.umn.edu), 
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca),
%          Benjamin VanderSluis (bvander@cs.umn.edu)
%
% Last revision: 2011-12-12
%
%%

function [all_linkage_cols, non_spec]  = filter_all_linkage_colonies_queryspecific_new(sgadata, linkagefile, all_querys, all_arrays, query_map, array_map, lfid)

    % Print the name and path of this script
    p = mfilename('fullpath');
    log_printf(lfid, '\nLinkage filter script:\n\t%s\n\n',p);
    
    % Load chromosomal coordinates and predefined linkage
    coord_fid = fopen('chrom_coordinates_111220.txt', 'r');
    coord_data = textscan(coord_fid, '%s%d%d%d', 'Delimiter', '\t', 'ReturnOnError', false);
    fclose(coord_fid);
    coord = struct();
    coord.strains = coord_data{1};
    coord.locations = [coord_data{3} coord_data{4}];
    coord.locations = sort(coord.locations, 2, 'ascend'); % enforce low to high assumption
    coord.locations = [coord_data{2} coord.locations];     % prepend chrom number
    
    linkage_fid = fopen(linkagefile, 'r');
    predef_lnkg_data = textscan(linkage_fid, '%s%d%d', 'Delimiter', '\t', 'ReturnOnError', false); 
    fclose(linkage_fid);
    predef_lnkg = struct();
    predef_lnkg.orf = predef_lnkg_data{1};
    predef_lnkg.coord = [predef_lnkg_data{2} predef_lnkg_data{3}];
    predef_lnkg.coord = sort(predef_lnkg.coord, 2, 'ascend'); % enforce low to high assumption

    % some stats for the log
    match_code_labels = {'query_strains_found_lnkg' 'query_orfs_found_lnkg' 'query_strains_found_coord' 'query_orfs_found_coord' 'querys_not_found'};
    match_code_counts = zeros(1,5);

    arrays_not_found_coord_exact = 0;
    arrays_not_found_coord_loose = 0;

    % ---------------------------------
    %% map array coordinates 
    array_coord = nan(length(all_arrays),3); % CHR START END
    for i = 1:length(all_arrays)
        array_orf = sgadata.orfnames{all_arrays(i)};

        ix = strmatch(array_orf, coord.strains, 'exact');

        if(isempty(ix))      % no exact match
            arrays_not_found_coord_exact = arrays_not_found_coord_exact + 1;
            ix = strmatch(strip_annotation(array_orf), coord.strains, 'exact');
        end

        if(isempty(ix))       % even the short version is missing
            arrays_not_found_coord_loose = arrays_not_found_coord_loose + 1;
            continue
        else
            array_coord(i,:) = coord.locations(ix,:);
        end
    end
    log_printf(lfid, 'Array coordinates mapped. \n');
    log_printf(lfid, '\texact match %d\n', length(all_arrays) - arrays_not_found_coord_exact);
    log_printf(lfid, '\tloose match %d\n', arrays_not_found_coord_exact - arrays_not_found_coord_loose);
    log_printf(lfid, '\tno match    %d\n', arrays_not_found_coord_loose);

    % ---------------------------------
    [ura3_linkage_arrays, ura3_match_code] = get_linked_arrays('YEL021W', predef_lnkg, coord, array_coord);
    [lyp1_linkage_arrays, lyp1_match_code] = get_linked_arrays('YNL268W', predef_lnkg, coord, array_coord);
    [can1_linkage_arrays, lyp1_match_code] = get_linked_arrays('YEL063C', predef_lnkg, coord, array_coord);
    lyp_can_linkage = unique([lyp1_linkage_arrays; can1_linkage_arrays]);

    all_linkage_cols_bool = boolean(zeros(size(sgadata.arrays))); % holds result, pre-allocated
    wild_type_id =  strmatch('undefined_sn4757', sgadata.orfnames);

    for i=1:length(ura3_linkage_arrays)
     for j=1:length(wild_type_id)
        all_linkage_cols_bool(intersect(array_map{all_arrays(ura3_linkage_arrays(i))}, query_map{wild_type_id(j)})) = true;
     end
    end
    for i=1:length(lyp_can_linkage)
        all_linkage_cols_bool(array_map{all_arrays(lyp_can_linkage(i))}) = true;
    end

    non_spec = find(all_linkage_cols_bool);

    % ---------------------------------
    %% determine specific linkage for each query 
    log_printf(lfid, ['Mapping query-specific linkage...\n|', blanks(50), '|\n|']);
    for i = 1:length(all_querys)
        strain_string = sgadata.orfnames{all_querys(i)};
        sgadata_this_query = sgadata.querys == all_querys(i);
        orf_list = split_orfs(strain_string);
        
        for j=1:length(orf_list)
            [this_orf_lnkg_arrays, match_code] = get_linked_arrays(orf_list{j},...
                predef_lnkg, coord, array_coord);
            match_code_counts(match_code) = match_code_counts(match_code)+1;

            for k = 1:length(this_orf_lnkg_arrays)
                all_linkage_cols_bool(intersect(array_map{all_arrays(this_orf_lnkg_arrays(k))}, query_map{all_querys(i)})) = true;
            end
        end

        % Print progress
        print_progress(lfid, length(all_querys), i);
    end

    log_printf(lfid, '|\nLinkage Query Process Report\n');
    for i=1:length(match_code_labels)
        log_printf(lfid, '\t%s\t%d\n', match_code_labels{i}, match_code_counts(i));
    end

    all_linkage_cols = find(all_linkage_cols_bool);
end
% ----------------------------------------------------------------------------------------------------------------------------------------
function [linked_arrays, match_code]  = get_linked_arrays(orf_string, predef_lnkg, coord, array_coord)
% returns indicies into all_arrays

    % Do we have predefined linkage for this orf?
    ix = strmatch(orf_string, predef_lnkg.orf, 'exact');

    % if not, do we have predified linkage for the unannotated version
    if(~isempty(ix))
        match_code = 1;
        left_win = predef_lnkg.coord(ix,1);
        right_win = predef_lnkg.coord(ix,2);
    else
        ix = strmatch(strip_annotation(orf_string), predef_lnkg.orf, 'exact');
    end

    % if still nothing, set the window from coordinates
    if(isempty(ix))
        [left_win, right_win, match_code] = lnkg_window_from_coord(orf_string, coord);
    else
        match_code = 2;
        left_win = predef_lnkg.coord(ix,1);
        right_win = predef_lnkg.coord(ix,2);
    end
    % At this point we have a linkage window defined for this orf
    % boolean, size of array (ASSUMPTION: upper case)
    this_orf_chrom = orf_string(2) - '@';
    % on this chrom, beginning, end (or both) within the window
    linked_arrays = find((array_coord(:,1) == this_orf_chrom) & ...
                          ((array_coord(:,2) > left_win & array_coord(:,2) < right_win) | ...
                           (array_coord(:,3) > left_win & array_coord(:,3) < right_win)));
end
% ----------------------------------------------------------------------------------------------------------------------------------------
function[left_win, right_win, match_code] = lnkg_window_from_coord(strain_string, coord)
    linkage_dist = 200e3;   % default linkage distance
    left_win = 1; % default values for error
    right_win = 1;
    ix = strmatch(strain_string, coord.strains, 'exact');
    if(~isempty(ix))
        match_code = 3;
    else
        ix = strmatch(strip_annotation(strain_string), coord.strains, 'exact');
        if(~isempty(ix))
            match_code = 4;
        else
            match_code = 5;
        end
    end

    if(~isempty(ix))
        left_win = max(coord.locations(ix,2) - linkage_dist, 1);
        right_win = coord.locations(ix,3) + linkage_dist;
    else
        return
    end
end
% ----------------------------------------------------------------------------------------------------------------------------------------
function[orf_list] = split_orfs(orf_string)
    ix = strfind(orf_string, '+');
    if(isempty(ix))
        orf_list = {orf_string};
    else
    	  orf_list = {orf_string(1:ix(1)-1) orf_string(ix(1)+1:end)};
    end
end 
