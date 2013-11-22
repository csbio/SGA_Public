function [sga] = sga_fill_square(strain_list, sga_sources)
%function [sga] = sga_fill_box(strain_list, sga_sources)
% given a list of querys, fill out a square with available epsilon and pvalue data
% queries are matched exactly, arrays are fuzzy
% enforces a reduce alleles on array axis
% asserts a reduce alleles on query axis


assert(length(strain_list) == length(unique(StripOrfs(strain_list))));


sga = struct();
sga.Cannon = struct();
sga.Cannon.Orf = StripOrfs(strain_list);
sga.Cannon.Common = OrfToCommon(StripOrfs(strain_list));
sga.Cannon.GENES = length(strain_list);
sga.Cannon.isArray = logical(ones(1,sga.Cannon.GENES));
sga.Cannon.isArray = logical(ones(sga.Cannon.GENES,1));

sga.eps = NaN(length(strain_list));
sga.pvl = NaN(length(strain_list));


% work backward and overwrite to give priority order

for i=length(sga_sources):-1:1
	source = reduce_alleles_array(sga_sources{i});

	% find the exact queries
	[Qint, qixa, qixb] = intersect(strain_list, source.Cannon.Orf(source.Cannon.isQuery));

	% fuzzy match the arrays
	[Aint, aixa, aixb] = intersect(StripOrfs(strain_list), ...
		StripOrfs(source.Cannon.Orf(source.Cannon.isArray)));


	S_EPS = source.eps(source.Cannon.isQuery, source.Cannon.isArray);
	S_PVL = source.pvl(source.Cannon.isQuery, source.Cannon.isArray);

	%fill them in ...???


