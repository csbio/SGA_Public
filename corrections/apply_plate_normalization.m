%%
% APPLY_PLATE_NORMALIZATION - corrects colony sizes for plate effects
%
% Inputs:
%   sgadata - structure containing all the colony size data
%   field - the name of the field to use as input
%   ignore_cols - indices of colonies to ignore in this correction
%   overall_med - reference colony size
%   plate_id_map - cell array of indices of colonies on each plate
%
% Outputs:
%   result - corrected colony sizes
%
% Authors: Chad Myers (cmyers@cs.umn.edu), Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
%
% Last revision: 2010-07-19
%
%%

function result = apply_plate_normalization(sgadata, field, ignore_cols, overall_med, plate_id_map)

    % Print the name and path of this script
    p = mfilename('fullpath');
    fprintf('\nPlate normalization script:\n\t%s\n\n',p);

    all_plates = unique(sgadata.plateids);
    result = sgadata.(field);
    vals  = sort(sgadata.(field)(~isnan(sgadata.(field))));

    if ~exist('overall_med')
        overall_med = nanmedian(vals(fix(.2*end):fix(.8*end)));
    end

    if ~exist('ignore_cols')
        ignore_cols = [];
    end

    savedata = sgadata.(field);
    sgadata.(field)(ignore_cols) = NaN;
    
    fprintf(['Plate normalization...\n|', blanks(50), '|\n|']);

    for i=1:length(all_plates)
        
        ind = plate_id_map{all_plates(i)};
        %ind = setdiff(ind,ignore_cols);  % not necessary if plate_id_map_ignorecols is passed.
   
        vals  = sort(sgadata.(field)(ind(~isnan(sgadata.(field)(ind)))));
        if length(vals) < 10
            continue;
        end

        plate_median = nanmedian(vals(fix(.2*end):fix(.8*end)));
        result(ind) = savedata(ind)*(overall_med/plate_median);
        
        % Print progress
        print_progress(length(all_plates), i);
        
    end
    fprintf('|\n');
