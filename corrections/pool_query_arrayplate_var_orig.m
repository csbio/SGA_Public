function [ query_arrplate_vars ] = pool_query_arrayplate_var_orig(sgadata, all_querys, query_map)
    query_arrplate_vars = [];   %zeros(length(all_querys),length(all_arrplates))+NaN;
    for i = 1:length(all_querys)
        ind = query_map{all_querys(i)};
        for j = 1:length(sgadata.all_arrayplateids)
            tmp_ind = find(sgadata.arrayplateids(ind) == sgadata.all_arrayplateids(j));
            currsets = unique(sgadata.setids(tmp_ind));
            for k = 1:length(currsets)
                tmp_ind2 = find(sgadata.setids(tmp_ind) == currsets(k));
                curr_ind = tmp_ind(tmp_ind2);
                query_arrplate_vars = [query_arrplate_vars; i, j, k,...
                   (1./(length(curr_ind)-length(unique(sgadata.arrays(curr_ind))))) * ...
                   nansum(sgadata.dm_normvar(curr_ind).*((sgadata.dm_num(curr_ind)-1)./sgadata.dm_num(curr_ind)))];
            end
        end
        
        % Print progress
        print_progress(lfid, length(all_querys), i);
    end
    log_printf(lfid, '|\n');
end