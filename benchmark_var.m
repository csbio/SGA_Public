

fprintf('alloc 1\n');
tic
A = pool_query_arrayplate_var_alloc(sgadata, all_querys, query_map, -11);
toc

fprintf('fixed 1\n');
tic
B = pool_query_arrayplate_var_fixed(sgadata, all_querys, query_map, -11);
toc

fprintf('alloc 2\n');
tic
A = pool_query_arrayplate_var_alloc(sgadata, all_querys, query_map, -11);
toc

fprintf('fixed 2\n');
tic
B = pool_query_arrayplate_var_fixed(sgadata, all_querys, query_map, -11);
toc
