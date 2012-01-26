function[residual, logresidual, arraymedian] = calculate_colony_residuals(sgadata, field, plate_id_map, lfid)

all_arrplates = unique(sgadata.arrayplateids);
width = 48;
height = 32;
% field = 'colsize_platenorm';

residual = sgadata.(field);
logresidual = sgadata.(field);
arraymedian = sgadata.(field);

array_means = zeros(length(all_arrplates),width*height)+NaN;

log_printf(lfid, ['Calculating colony residuals...\n|' blanks(50) '|\n|']);

for i = 1:length(all_arrplates)
    currplates = unique(sgadata.plateids(sgadata.arrayplateids == all_arrplates(i)));
    t=[];   % stores the list of colony sizes for this arrayplate
    save_inds = struct;

    for j = 1:length(currplates)

        % Get a matrix of colony sizes (d) and colony indices (d_map) for the plate
        ind = plate_id_map{currplates(j)};
        d = zeros(32,48);
        d(:,[1,2,47,48]) = NaN;
        d([1,2,31,32],:) = NaN;
 
        d_map = zeros(32,48)+1; % Default reference to 1st index

        iii = sub2ind([32,48], sgadata.rows(ind), sgadata.cols(ind));
        d_map(iii) = ind;
        d(iii) = sgadata.(field)(ind);

        t = [t,d(:)];
        save_inds(j).ind_mat = d_map(:);
        save_inds(j).dat_mat = d(:);

    end

    array_means(i,:) = nanmedian(t,2);
    array_means(i,array_means(i,:) <= 0)=NaN;

    for j=1:length(currplates)

        t = sgadata.(field)(save_inds(j).ind_mat);
        t(t <= 0)=NaN;
        logresidual(save_inds(j).ind_mat)=log(t./array_means(i,:)');
        arraymedian(save_inds(j).ind_mat)=array_means(i,:)';
        residual(save_inds(j).ind_mat)=t-array_means(i,:)';

    end

    % Print progress
    print_progress(lfid, length(all_arrplates), i);

end

log_printf(lfid, '|\n'); 


%{
% Old method
% Filter very large colonies
ind = find(sgadata.colsize_platenorm >= 1.5*default_median_colsize & ...
    sgadata.rows > 2 & sgadata.rows < 31 & ...
    sgadata.cols > 2 & sgadata.cols < 47);

all_spots = unique(sgadata.spots(ind));
num_big_colonies_per_spot = histc(sgadata.spots(ind), all_spots);
spots_to_remove = all_spots(num_big_colonies_per_spot >= 3);

ii = find(ismember(sgadata.spots, spots_to_remove));
sgadata.colsize_platenorm(ii) = NaN;
%}
% 
