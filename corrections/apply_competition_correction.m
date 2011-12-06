%%
% APPLY_COMPETITION_CORRECTION - corrects colony sizes for competition effects
%
% Inputs:
%   sgadata - structure containing all the colony size data
%   field - the name of the field to be used as input
%   ignore_cols - indices of colonies to ignore in this correction
%   plate_id_map - cell array of indices of colonies on each plate
%
% Outputs:
%   newdata - corrected colony sizes
%
% Authors: Chad Myers (cmyers@cs.umn.edu), Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
%
% Last revision: 2010-07-19
%
%%

function newdata = apply_competition_correction(sgadata,field,ignore_cols,plate_id_map,lfid)

    % Print the name and path of this script
    p = mfilename('fullpath');
    log_printf(lfid, '\nCompetition correction script:\n\t%s\n',p);    
    
    % Get indices of neighboring colonies
    colony_neighbor_inds = get_colony_neighbor_indices_list(lfid);
    
    % List of unique plateids
    all_plates = unique(sgadata.plateids);
    
    % List of colony sizes for each colony's neighbors
    sgadata.neighbor_cols = zeros(length(sgadata.colsize),3);   

    log_printf(lfid, ['Mapping neighboring colonies...\n|', blanks(50), '|\n|']);
    for i=1:length(all_plates)
        
        % List of sgadata indices to this plate colonies
        ind = plate_id_map{all_plates(i)};
        
        % Array of sgadata indices to this plate colonies
        plate_map = zeros(32,48)+1;  %default reference to 1st index
        
        % Populate the array
        iii = sub2ind([32,48], sgadata.rows(ind),sgadata.cols(ind));
        plate_map(iii) = ind;
       
        % Build the list of sgadata indices to this plate colony neighbors
        nhbr_inds = zeros(size(colony_neighbor_inds));
        nhbr_inds(:) = plate_map(colony_neighbor_inds(:));
        
        % Store colony neighbors sizes to sgadata
        sgadata.neighbor_cols(plate_map(:),:) = sgadata.(field)(nhbr_inds);
        
        % Print progress
        print_progress(lfid, length(all_plates),i);
        
    end 
    log_printf(lfid, '|\n');
    
    sgadata.neighbor_cols(sgadata.neighbor_cols == sgadata.(field)(1)) = NaN;
    
    % Mask the borders
    inds = find(sgadata.rows < 3 | sgadata.rows > 30 | sgadata.cols < 3 | sgadata.cols > 46);
    sgadata.neighbor_cols(inds,:) = 0;

    
    %% Do quantile-quantile normalization on colony size distribution binned by neighbor size
    
    keep_ind = find(nansum(sgadata.neighbor_cols,2) ~= 0);
    [vals,sort_ind] = sort(min(sgadata.neighbor_cols(keep_ind,:),[],2));
    col_cutoffs = [-inf;vals(fix(end*[.1:.1:1]))];

    % Get quantile of background distribution (middle 60-80% neighbor sizes)
    tmp = sgadata.(field)(keep_ind(sort_ind(fix(end*.6):fix(end*.8))));
    tmp(isnan(tmp))=[];

    comp_corr_global=zeros(length(sgadata.colsize),1);
    newdata = sgadata.(field);
    
    log_printf(lfid, ['\nCompetition correction...\n|' blanks(50) '|\n|']);
    for i = 2:(fix(length(col_cutoffs)/2)+2) % correct 0-70% bins each independently
        
        ind2 = find(min(sgadata.neighbor_cols,[],2) >= col_cutoffs(i-1) & ...
            min(sgadata.neighbor_cols,[],2) < col_cutoffs(i) & ...
            nansum(sgadata.neighbor_cols,2) ~= 0 & ...
            ~ismember([1:length(sgadata.(field))]',ignore_cols));
        
        
%         if length(ind2) < 100
%             continue;
%         end

        currdata = sgadata.(field)(ind2);   
        
        % Do within-group smoothing first
        [v,ind] = sort(min(sgadata.neighbor_cols(ind2,:),[],2),'ascend');
        lowsmoothed = smooth(min(sgadata.neighbor_cols(ind2(ind),:),[],2),double(currdata(ind)),500,'lowess');

        comp_corr_global(ind2(ind)) = lowsmoothed;
        
        % Now, do quantile normalization on upper half
        data2 = quantile_normalization(currdata-comp_corr_global(ind2),tmp);
        
        % Get complete correction
        newdata(ind2) = data2;
        
        % Print progress
        print_progress(lfid, length(2:(fix(length(col_cutoffs)/2)+2)), i);
        
    end
    log_printf(lfid, '|\n');
    
    newdata = sgadata.rowcolcorr_colsize - (sgadata.residual_spatialnorm - newdata);
    
