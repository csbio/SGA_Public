
% These are the defaults shared by all jobs
param = struct();
param.skip_wt_remove = true;
param.skip_linkage_detection = true;
param.wild_type = 'URA3control_sn4757';
param.border_strain_orf = 'YMR271C_dma3865';
param.linkagefile = 'refdata/linkage-est_sga_merged_131122.txt';
param.coord_file = 'refdata/chrom_coordinates_111220.txt';
param.removearraylist = 'refdata/bad_array_strains_140526.csv';

save rawdata/chi/141027/job_defaults.mat param


% examples
% inputfile = 'rawdata/chi/141027/collab_dump_141027_CHI_Michael_T30_day1_SC_setA.txt';
% outputfile= 'scored/chi/141027/scored_chi_t30_day1_141027_scored_141001';
% smfitnessfile = 'refdata/WTonly_chi_ts_t30_120524_SMF_120627_arrayAVE.txt';

replicates = {'', 'replicate/'};
temps = {'T30', 'T38'};
days = {'day1', 'day2'};
medias = {'SC', 'YPED'};
sets = {'A', 'B', 'C'};
counter = 1;

for r=1:2
for t=1:2
for d=1:2
for m=1:2
for s=1:3

    % reset param
    load rawdata/chi/141027/job_defaults.mat param

    % set file specific parameters
    param.inputfile = sprintf('rawdata/chi/141027/%scollab_dump_141027_CHI_Michael_%s_%s_%s_set%s.txt',...
        replicates{r}, temps{t}, days{d}, medias{m}, sets{s});
    fprintf('******************** Starting Job %d / 48:\n%s\n\n', counter, param.inputfile);

    param.outputfile= sprintf('scored/chi/141027/%sscored_chi_1plate_141027_%s_%s_%s_set%s_scored_141112',...
        replicates{r}, temps{t}, days{d}, medias{m}, sets{s});

    param.smfitnessfile = sprintf('refdata/WTonly_chi_ts_%s_120524_SMF_120627_arrayAVE.txt', lower(temps{t}));

    % run the job
    compute_sgascore_function(param);
    counter = counter+1;
end
end
end
end
end
