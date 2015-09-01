    f=1;
    path_parts = split('/', inputfile);
    datasets(f).name = path_parts{end};
    datasets(f).path = inputfile;
    
    [datasets(f).queries, datasets(f).arrays, datasets(f).scores] = read_matrix_file([inputfile, '_col2.txt'],2,2);
    [datasets(f).queries, datasets(f).arrays, datasets(f).scores_sd] = read_matrix_file([inputfile, '_col3.txt'],2,2);
    [datasets(f).queries, datasets(f).arrays, datasets(f).scores_pvalue] = read_matrix_file([inputfile, '_col4.txt'],2,2);
   
    [datasets(f).queries, datasets(f).arrays, datasets(f).query_smf_matrix] = read_matrix_file([inputfile, '_col5.txt'],2,2);
    [datasets(f).queries, datasets(f).arrays, datasets(f).query_smf_matrix_sd] = read_matrix_file([inputfile, '_col6.txt'],2,2);
    [datasets(f).queries, datasets(f).arrays, datasets(f).array_smf_matrix] = read_matrix_file([inputfile, '_col7.txt'],2,2);
    [datasets(f).queries, datasets(f).arrays, datasets(f).array_smf_matrix_sd] = read_matrix_file([inputfile, '_col8.txt'],2,2);
    
    
  
    
    [datasets(f).queries, datasets(f).arrays, datasets(f).dm_obs] = read_matrix_file([inputfile, '_col10.txt'],2,2);
    [datasets(f).queries, datasets(f).arrays, datasets(f).dm_obs_sd] = read_matrix_file([inputfile, '_col11.txt'],2,2);
    
    [datasets(f).queries, datasets(f).arrays, datasets(f).dm_exp] = read_matrix_file([inputfile, '_col9.txt'],2,2);
    
    datasets(f).scores_eps = datasets(f).dm_obs - datasets(f).dm_exp;
    datasets(f).scores_eps_sd = datasets(f).dm_obs_sd;
    
    
    
    [datasets(f).queries, datasets(f).arrays, datasets(f).bckg_mean] = read_matrix_file([inputfile, '_col12.txt'],2,2);
    [datasets(f).queries, datasets(f).arrays, datasets(f).bckg_std] = read_matrix_file([inputfile,'_col13.txt'],2,2);
    
    inds = find(datasets(f).scores == 0 | isnan(datasets(f).scores) | datasets(f).scores_sd == 0 | isnan(datasets(f).scores_sd) | isnan(datasets(f).scores_pvalue));
    
    datasets(f).scores(inds) = NaN;
    datasets(f).scores_sd(inds) = NaN;
    datasets(f).scores_pvalue(inds) = NaN;
    datasets(f).dm_obs(inds) = NaN;
    datasets(f).dm_obs_sd(inds) = NaN;
    datasets(f).dm_exp(inds) = NaN;
    datasets(f).scores_eps(inds) = NaN;
    datasets(f).scores_eps_sd(inds) = NaN;
    



    [t,ix1] = sort_by_chr(datasets(f).queries);
    [t,ix2] = sort_by_chr(datasets(f).arrays);




    datasets(f).queries = datasets(f).queries(ix1);
    datasets(f).arrays = datasets(f).arrays(ix2);
    
    flds = fieldnames(datasets(f));
    for i = 5 : length(flds)
        datasets(f).(flds{i}) = datasets(f).(flds{i})(ix1,ix2);
    end
    
    
    
    
    
    
%     pvals = sqrt(normcdf(-abs(fg_120118.scores_corr./fg_120118.scores_sd)) .* normcdf(-abs(log((fg_120118.bckg_mean + fg_120118.scores_corr)./fg_120118.bckg_mean) ./ log(fg_120118.bckg_std) )));

