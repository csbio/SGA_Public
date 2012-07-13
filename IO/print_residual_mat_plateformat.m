function [big_map, labels_y, labels_x] = print_residual_mat_plateformat_rewrite(sgadata,field,filename, plate_id_map, center_arrays)
%function [big_map, labels_y, labels_x] = print_residual_mat_plateformat_rewrite(...
%     sgadata, field, filename, plate_id_map, center_arrays[T/F])
% Benjamin VanderSluis (bvander@cs.umn.edu)
% April 10, 2012
total_time = tic();
all_arrplates = unique(sgadata.arrayplateids);
width = 48;
height = 32;
extra_border = 0; % add this much around each plate  (1 = 1 col between plates)
border_mask = 2;  % mask out this much for each plate (1 = 1 col from each plate, 2 between)

if(~exist('center_arrays', 'var'))
	center_arrays = false;
end


% build a list of queryx set pairs, number of plate rows
plate_rows = unique([sgadata.querys sgadata.setids], 'rows');
plate_centers = cell(1,length(all_arrplates));
	plate_centers(:) = {zeros(height, width)};
plate_counts  = cell(1, length(all_arrplates));
	plate_counts(:) = {zeros(height, width)};

% a place to hold our result
big_size_y = length(plate_rows)*(height+extra_border)+extra_border;
big_size_x = length(all_arrplates)*(width+extra_border)+extra_border;

big_map = nan(big_size_y, big_size_x);
	labels_y = cell(big_size_y,1);
	labels_x = cell(big_size_x,1);
	labels_y(:) = {'brdr'};
	labels_x(:) = {'brdr'};

% precalculate array plate (col) matches
plate_columns = boolean(zeros(length(sgadata.arrayplateids), length(all_arrplates)));
for i=1:length(all_arrplates)
	plate_columns(:,i) = sgadata.arrayplateids == all_arrplates(i);
end

% main loop
for r=1:length(plate_rows)
	this_row_bool = (sgadata.querys == plate_rows(r,1)) & (sgadata.setids == plate_rows(r,2));
	for c = 1:length(all_arrplates)
		
		[big_r, big_c] = calculate_global_index(r, c, extra_border, height, width);

		% update the labels % before errors
		for i=1:height
			labels_y{big_r+i-1} = ['Query: ',sgadata.orfnames{plate_rows(r,1)},' Set: ', num2str(plate_rows(r,2)), ' Row: ',num2str(i)];
		end

		% find THIS plate
		plate_ix = find(this_row_bool & plate_columns(:,c));
		if(length(plate_ix) > 0) % skip empty plates
			if(length(plate_ix) > height*width)
				fprintf('collision Query: %s\tPlate: %d\tSet %d\n', sgadata.orfnames{plate_rows(r,1)}, c, plate_rows(r,2));
				continue;
			end

			plate_data = nan(height, width);

			plate_row = sgadata.rows(plate_ix);
			plate_col = sgadata.cols(plate_ix);
			plate_val = sgadata.(field)(plate_ix);

			plate_ind = sub2ind([height, width], plate_row, plate_col);
			plate_data(plate_ind) = plate_val(:);

			plate_counts{c} = plate_counts{c}+~isnan(plate_data);
			plate_centers{c}(~isnan(plate_data)) = plate_centers{c}(~isnan(plate_data)) + plate_data(~isnan(plate_data));
			
			if(border_mask > 0)
				plate_data(1:border_mask,:) = NaN;
				plate_data(end-border_mask+1:end,:) = NaN;
				plate_data(:,1:border_mask) = NaN;
				plate_data(:,end-border_mask+1:end) = NaN;
			end

			% put this plate in its place on the whole map
			big_map(big_r:big_r+height-1, big_c:big_c+width-1) = plate_data;
		end

	end
end

if(center_arrays)
	% Lets add a fake query representing the centers
	big_map = [big_map; nan(height+extra_border, (width+extra_border)*length(all_arrplates)+extra_border)];
	labels_y(end+1: size(big_map,1)) = {'center'};

	for c=1:length(all_arrplates)
		plate_centers{c} = plate_centers{c} ./ plate_counts{c};
		for r=1:length(plate_rows)
			[big_r, big_c] = calculate_global_index(r, c, extra_border, height, width);
			big_map(big_r:big_r+height-1, big_c:big_c+width-1) = big_map(big_r:big_r+height-1, big_c:big_c+width-1) - plate_centers{c};
		end
	end

	for c=1:length(all_arrplates)
		[big_r, big_c] = calculate_global_index(r+1, c, extra_border, height, width);
		big_map(big_r:big_r+height-1, big_c:big_c+width-1) = plate_centers{c};
	end
	
end

		
% y axis labels
ptr = 1;
for j=1:length(all_arrplates),
	ptr = ptr+extra_border;
    for i=1:width,
        labels_x{ptr} = ['Array plate ',num2str(j),'; Column ',num2str(i)];
        ptr=ptr+1;
    end
end

fprintf('data tiling complete, begin writing outputfile\n');
print_pcl_file(big_map, labels_y, labels_y, labels_x, filename);
tt = toc(total_time);

fprintf('total time elapsed: %.2fs\n', tt);

end % MAIN

function[global_r, global_c] = calculate_global_index(query, arrayplate, extra_border, height, width)
	% returns the global origin for a single plate (upper left corner)
	% ROWS = num_queries*(height + extra_border) + extra_border (on the bottom)
	% COLS = num_array_plates *(width + extra_border + extra_border (on the end)
	global_r = (query-1)     *(height + extra_border) + extra_border + 1;
	global_c = (arrayplate-1)*(width  + extra_border) + extra_border + 1;
end
