function [ query_arrplate_vars ] = pool_query_arrayplate_var(sgadata, all_querys, query_map, lfid)
    query_arrplate_vars = [];
    for i = 1:length(all_querys)
        % colonies from _this_ query
        % q_ix -> sgadata()
        q_ix = query_map{all_querys(i)}; 

        for j = 1:length(sgadata.all_arrayplateids)
            % qa_ix finds colonies for _this_ query on _this_ arrayplate
            % qa_ix -> sgadata(q_ix);
            qa_ix = find(sgadata.arrayplateids(q_ix) == sgadata.all_arrayplateids(j));
            qa_sets = unique(sgadata.setids(q_ix(qa_ix)));

            for k = 1:length(qa_sets)
                tmp_ind2 = find(sgadata.setids(q_ix(qa_ix)) == qa_sets(k));
                qas_ix = q_ix(qa_ix(tmp_ind2));

                num_cols = length(qas_ix); % total colonies, this query, this arrayplate, this set
                num_arrays = length(unique(sgadata.arrays(qas_ix)));
                df = 1 / (num_cols - num_arrays); % scalar, degrees of freedom?
                xxx = nansum(sgadata.dm_normvar(qas_ix) ...
                         .*  (sgadata.dm_num(qas_ix)-1) ...
                         ./  sgadata.dm_num(qas_ix));

                var_row = [i, j, k, df * xxx];
                query_arrplate_vars = [query_arrplate_vars; var_row];
            end
        end
        
        % Print progress
        print_progress(lfid, length(all_querys), i);
    end
    log_printf(lfid, '|\n');
end
