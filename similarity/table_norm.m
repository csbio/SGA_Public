function [ data3_norm, table] = table_norm(data1, data2, data3)
%function [ data3_norm, table, data2_norm ] = table_norm(data1, data2, data3)
% normaliza data according to a table 
% data1 and data2 are matched observations
% 		data2 gets quantile_normed to data1
% 		data3 gets adjusted based on data1:data2 reference


	assert(length(data1) == length(data2), 'First two arguments must be matched observations.');

	% Deal with NaNs
	nn = ~isnan(data1) & ~isnan(data2);
	data2_norm = NaN(size(data2));
	data2_norm(nn) = quantile_normalization(data2(nn), data1(nn));



	% translation table maps old data2 values to new ones
	% sort ix is the same pre and post...
	table = [data2(nn) data2_norm(nn)];
	[~, ix] = sort(data2(nn), 'ascend');
	table = table(ix,:);

	% remove redundant key entries
	% keep the median of the values
	[t1, head] = unique(table(:,1), 'first');
	[~, tail] = unique(table(:,1), 'last');
	t2 = NaN(size(t1));

	for i=1:length(t2)
		t2(i) = median(table(head(i):tail(i), 2));
	end

	table = [t1 t2];

	% table index = how many table values are < data
	data3_norm = NaN(size(data3));
	data3_nn = data3(~isnan(data3));

	% so we'll do one table row at a time
	% fast but all at once takes too much memory
	data3_tableix = zeros(size(data3_nn));
% 	progress = 0;
	N = length(t1);
	fprintf('\n\n\n\n');
	for i=1:N
		data3_tableix = data3_tableix + double(data3_nn > t1(i));
% 		progress = printprogress(i, N, progress);
	end

	data3_tableix(data3_tableix == 0) = 1; % more extreme than we saw or NaN

	data3_norm(~isnan(data3)) = t2(data3_tableix);
end

