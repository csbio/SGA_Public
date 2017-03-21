%%
% QUANTILE_NORMALIZATION - normalizes the distribution to a reference
%
% Inputs:
%   data - data to be normalized
%   refdist - reference distribution
%
% Outputs:
%   norm_dist - normalized distribution
%
% Authors: Chad Myers (cmyers@cs.umn.edu), Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
%
% Last revision: 2010-07-19
%
%%

function norm_dist = quantile_normalization(data,refdist)

ind = find(~isnan(data));

quantiles = 1/length(ind):1/length(ind):1;
refquant = quantile(refdist,quantiles);

[vals,sort_ind]=sort(data(ind),'ascend');
norm_dist = data;
norm_dist(ind(sort_ind))=refquant;
