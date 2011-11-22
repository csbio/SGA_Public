%%
% APPLY_JACKKNIFE_CORRECTION - removes colonies which contribute hugely to variance
%
% Inputs:
%	 sgadata - structure containing all the colony size data
%	 field - the name of the field to be used as input
%	 plate_id_map - cell array of indices of colonies on each plate
%
% Outputs:
%	 newdata - corrected colony sizes
%
% Authors: Chad Myers (cmyers@cs.umn.edu)
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
%          Benjamin VanderSluis (bvander@cs.umn.edu)
%
% Last revision: 2011-11_21
%
%%

function result = apply_jackknife_correction(sgadata,field,border_strain_orf,query_map,plate_id_map)

	% Print the name and path of this script
	p = mfilename('fullpath');
	fprintf('\nCompetition correction script:\n\t%s\n',p);		

	all_querys = unique(sgadata.querys);
	result = sgadata.(field);

	his3_ind = strmatch(border_strain_orf,sgadata.orfnames,'exact');
	assert(~isempty(his3_ind), 'is %s the correct border strain?\n', border_strain_orf);

	fprintf(['Running the hold-one-out filter...\n|' blanks(50) '|\n|']);
	for i = 1:length(all_querys)
		
		ind1 = query_map{all_querys(i)};
		curr_plates = unique(sgadata.plateids(ind1));
		
		for k = 1:length(curr_plates)
				
			ind = plate_id_map{curr_plates(k)};
			max_cols = length(find(sgadata.arrays(ind) == his3_ind));

			curr_arrays = unique(sgadata.arrays(ind));
			curr_mat = zeros(length(curr_arrays),max_cols)+NaN;
			ind_mat = zeros(length(curr_arrays),max_cols)+NaN;
	
			for j = 1:length(curr_arrays)	 
				ind2 = find(sgadata.arrays(ind)==curr_arrays(j));
				curr_mat(j,1:length(ind2)) = sgadata.(field)(ind(ind2));
				ind_mat(j,1:length(ind2)) = ind(ind2);

				if curr_arrays(j) == his3_ind
					continue;
				end
		
				ind3 = find(~isnan(curr_mat(j,:)));
				vals = jackknife(@nanstd,curr_mat(j,ind3));		 
						
				tot_dev = nanvar(curr_mat(j,ind3))*(length(vals)-1);
				jackknife_dev = (vals.^2).*(length(vals)-2);
						
				t = find(tot_dev - jackknife_dev > 0.9*tot_dev); % find colonies that contribute more than 90% of total variance.
						
				if length(t) <= 0.25 * length(vals)						
					result(ind_mat(j,ind3(t))) = NaN;
				end

			end
				
		end
	 
		% Print progress
		print_progress(length(all_querys), i);

	end
	fprintf('|\n');

end


%{
     x = curr_mat(j,ind3);
     vars = (x-nanmean(x)).^2;
     
            t = find(vars./sum(vars) > 0.9); % find colonies that contribute more than 90% of total variance.
%}
