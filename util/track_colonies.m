function[results, cell_results] = track_colonies(query_list, array_list, sgadata, all_querys, all_arrays, linkage_cols, bad_array_cols)
%function[results, cell_results] = track_colonies(query_list, array_list, sgadata, all_querys, all_arrays, linkage_cols, bad_array_cols)
% designed to track a %SMALL% number of colonies through the pipleine

q_ix = [];
for i=1:length(query_list)
	q_ix = [q_ix; strmatch(query_list{i}, sgadata.orfnames, 'exact')];
end

a_ix = [];
for i=1:length(array_list)
	a_ix = [a_ix; strmatch(array_list{i}, sgadata.orfnames)];
	%a_ix = [a_ix; strmatch([array_list{i} '_dma'], sgadata.orfnames)];
end

%% code for all pairs

colonies_to_watch = find(ismember(sgadata.arrays, a_ix) & ismember(sgadata.querys, q_ix));


%% code for matched pairs
%assert(length(q_ix) == length(a_ix));
%colonies_to_watch = [];
%for i=1:length(q_ix)
%	colonies_to_watch = [colonies_to_watch; find(sgadata.arrays == a_ix(i) & sgadata.querys == q_ix(i))];
%end

% not all of these are "colony" scores, some are repmat(array_scores e.g.) 
% but they are in "chronological" order
fields_to_track = {...
	'colsize',...
	'colsize_platenorm',...
	'residual',...
	'logresidual',...
	'arraymedian',...
	'spatialnorm_colsize',...
	'rowcolcorr_colsize',...
	'residual_spatialnorm',...
	'finalplatecorr_colsize',...
	'filt_jackknife',...
	'filt_colsize',...
	'batchnorm_colsize',...
	'batchnorm_colsize_nonegs',...
	'arraymean_corrected'};

	%'compcorr_colsize',...


results = zeros(length(colonies_to_watch), length(fields_to_track));
for f=1:length(fields_to_track)
	results(:,f) = sgadata.(fields_to_track{f})(colonies_to_watch);
end

%fields_to_track = [fields_to_track, 'linkage_cols', 'bad_array_cols'];
%results = [results zeros(size(results,1),2)];
%results(:,f+1) = ismember(colonies_to_watch, linkage_cols);
%results(:,f+2) = ismember(colonies_to_watch, bad_array_cols);


cell_results = cell(size(results,1)+1, size(results,2)+2);
cell_results(2:end, 3:end) = num2cell(results);
cell_results(2:end, 1) = sgadata.orfnames(sgadata.querys(colonies_to_watch));
cell_results(2:end, 2) = sgadata.orfnames(sgadata.arrays(colonies_to_watch));
cell_results(1,3:end) = fields_to_track;
cell_results(1,[1,2]) = {'Query', 'Array'};



