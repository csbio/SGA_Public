function [ query_arrplate_vars ] = pool_query_arrayplate_var(sgadata, all_querys, query_map, lfid)
    query_arrplate_vars = [];   %zeros(length(all_querys),length(all_arrplates))+NaN;
    for i = 1:length(all_querys)
        % colonies from _this_ query
        % q_ix -> sgadata()
        q_ix = query_map{all_querys(i)}; 

        for j = 1:length(sgadata.all_arrayplateids)
            % qa_ix finds colonies for _this_ query on _this_ arrayplate
            % qa_ix -> sgadata(q_ix);
            qa_ix = find(sgadata.arrayplateids(q_ix) == sgadata.all_arrayplateids(j));

            % Here we make the first error:
            % qa_sets is the list of set (replicate) ids found for 
            % _this_ query on _this_ arrayplate
            % qa_ix -> sgadata(q_ix), NOT sgadata()
            qa_sets = unique(sgadata.setids(qa_ix));
            % as a result we are not even iterating over the right sets
            % so all of the values are wrong, and MOST of the values are 0

            % correct would be
            % qa_sets = unique(sgadata.setids(q_ix(qa_ix)));

            for k = 1:length(qa_sets)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % qas_ix finds colonies for _this_ set (and this q/ap)
                % qas_ix -> sgadata()
                % the problem is sgadata.setids(qa_ix) is nonsense 
                % because qa_ix -> sgadata(q_ix)
                %    NOT  qa_ix -> sgadata()
                tmp_ind2 = find(sgadata.setids(qa_ix) == qa_sets(k));
                qas_ix = qa_ix(tmp_ind2);
                % the damage is done, and the second line cannot unscramble it

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % this would be one correct way
                % tmp_ind2 = find(sgadata.setids(q_ix(qa_ix)) == qa_sets(k));
                % qas_ix = q_ix(qa_ix(tmp_ind2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
