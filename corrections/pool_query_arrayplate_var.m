function query_arrplate_vars = pool_query_arrayplate_var(sgadata, all_querys, query_map, lfid)

    % this is a generous guess as to how many rows we'll need 
    % enough for 10 reps of each query, but still probably < 1M
    N = length(all_querys) * length(sgadata.all_arrayplateids) * 10;

    query_arrplate_vars = NaN(N,4);
    ptr = 1; % keep track of where we are and trim the tail

    for i = 1:length(all_querys)
        % colonies from _this_ query
        % q_ix -> sgadata()
        q_ix = query_map{all_querys(i)}; 

        for j = 1:length(sgadata.all_arrayplateids)
            % qa_ix finds colonies for _this_ query on _this_ arrayplate
            % qa_ix -> sgadata(q_ix);
            qa_ix = find(sgadata.arrayplateids(q_ix) == ...
                         sgadata.all_arrayplateids(j));

            qa_sets = unique(sgadata.setids(q_ix(qa_ix)));

            for k = 1:length(qa_sets)
                tmp_ind2 = find(sgadata.setids(q_ix(qa_ix)) == qa_sets(k));
                qas_ix = q_ix(qa_ix(tmp_ind2));

                num_cols = length(qas_ix); % colonies for this query-arrayplate-set
                num_arrays = length(unique(sgadata.arrays(qas_ix)));
                df = 1 / (num_cols - num_arrays); % scalar, degrees of freedom
                xxx = nansum(sgadata.dm_normvar(qas_ix) ...
                         .*  (sgadata.dm_num(qas_ix)-1) ...
                         ./  sgadata.dm_num(qas_ix));

                if ptr > N
                    % we need another chunk of space, if you see this all the time
                    % our guess was bad. Bump up N at the top.
                    log_printf(' -- query_arrplate_vars space reallocation --\n');
                    query_arrplate_vars = [query_arrplate_vars; NaN(N,4)]
                    N = 2*N;
                end

                query_arrplate_vars(ptr,:) = [i, j, k, df * xxx];
                ptr = ptr+1;
            end
        end
        
        print_progress(lfid, length(all_querys), i);
    end
    query_arrplate_vars = query_arrplate_vars(1:ptr-1,:); % trim exess nans
    log_printf(lfid, '|\n');
end
