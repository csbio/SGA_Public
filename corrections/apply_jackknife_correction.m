%%
% APPLY_JACKKNIFE_CORRECTION - removes colonies which contribute hugely to variance
%
% Inputs:
%	 sgadata - structure containing all the colony size data
%	 field - the name of the field to be used as input
%   border_strain_orf - name of the border strain to be skipped
%	 plate_id_map - cell array of indices of colonies on each plate
%
% Outputs:
%	 newdata - corrected colony sizes
%
% Authors: Chad Myers (cmyers@cs.umn.edu)
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
%          Benjamin VanderSluis (bvander@cs.umn.edu)
%
% Last revision: 2011-11-29
%
%%

function result = apply_jackknife_correction(sgadata,field,border_strain_orf,query_map,plate_id_map,lfid)

	% Print the name and path of this script
	p = mfilename('fullpath');
	log_printf(lfid, '\nJackknife variance correction script:\n\t%s\n\n',p);		

	all_querys = unique(sgadata.querys);
	result = sgadata.(field);

	border_id = strmatch(border_strain_orf,sgadata.orfnames, 'exact');
	assert(length(border_id) == 1, 'is %s the correct border strain?\n', border_strain_orf);

	log_printf(lfid, ['Running the hold-one-out filter...\n|' blanks(50) '|\n|']);
	for i = 1:length(all_querys)
		
		ind1 = query_map{all_querys(i)};
		curr_plates = unique(sgadata.plateids(ind1));
		
		for k = 1:length(curr_plates)
				
			ind = plate_id_map{curr_plates(k)};
			max_cols = sum(sgadata.arrays(ind) == border_id);

			curr_arrays = unique(sgadata.arrays(ind));
			curr_mat = zeros(length(curr_arrays),max_cols)+NaN;
			ind_mat = zeros(length(curr_arrays),max_cols)+NaN;
	
			for j = 1:length(curr_arrays)	 
				ind2 = find(sgadata.arrays(ind)==curr_arrays(j));
				curr_mat(j,1:length(ind2)) = sgadata.(field)(ind(ind2));
				ind_mat(j,1:length(ind2)) = ind(ind2);

				if curr_arrays(j) == border_id
					continue;
				end
		


				ind3 = find(~isnan(curr_mat(j,:)));

				tot_var = var(curr_mat(j,ind3))*(length(ind3)-1);
				jack_std = jackknife(@std,curr_mat(j,ind3));		 
				jackknife_dev = (jack_std.^2).*(length(jack_std)-2);
				t = find(jackknife_dev < 0.1*tot_var); % find colonies that contribute more than 90% of total variance

     			%x = curr_mat(j,ind3);
     			%sqrd_devs = (x-mean(x)).^2;
            %t = find(sqrd_devs ./ sum(sqrd_devs) > 0.75); 
						
				if length(t) <= 0.25 * length(ind3)						
					result(ind_mat(j,ind3(t))) = NaN;
				end

			end
				
		end
	 
		% Print progress
		print_progress(length(all_querys), i);

	end
	log_printf(lfid, '|\n');

end

