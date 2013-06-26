function[eps] = quantile_normalize_from_table(eps, table)
%function[eps] = quantile_normalize_from_table(eps, table)
% uses a table of cannonical values from a previous, matched distribution
% quantile normalization, to adjust the distribution of values in eps
% Benjamin VanderSluis (bvander@cs.umn.edu)
% 130603

	data = eps(~isnan(eps));
	data_tableix = zeros(size(data));

	for i=1:size(table,1)
		% count how many values in col 1 are less than each data point
		data_tableix = data_tableix + double(data > table(i,1));
	end

	data_tableix(data_tableix == 0) = 1; 
	eps(~isnan(eps)) = table(data_tableix,2);
end
