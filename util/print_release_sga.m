function[] = print_release_sga(sga, outputfile, Experiment);
%function[] = print_release_sga(sga, outputfile, Experiment);
% prints a release file from a struct. Include an EXP string (e.g. TS26)
% Can print EXP from struct if available
	%{
	----------------------------------------------------------
	TYPE RELEASE
	9 col
	----------------------------------------------------------
	1  Query ORF 
	2  Query gene name
	3  Array ORF 
	4  Array gene name
	5  Double mutant fitness
	6  Genetic interaction score (eps)
	7  Standard deviation
	8  P-value
	9  Experiment
	%}

	if ~isfield(sga, 'src')
		sga.src = ones(sga.Cannon.GENES);
		sga.sources = {Experiment};
	end

	fid = fopen(outputfile, 'w');

	header = {'Query ORF', 'Query gene name', 'Array ORF', 'Array gene name', ...
	          'Double mutant fitness', 'Epsilon', 'Standard deviation', ...
	          'P-value', 'Experiment'};
	for i=1:length(header)-1
		fprintf(fid, '%s\t', header{i});
	end
	fprintf(fid, '%s\n', header{end});

	[r, c] = find(~isnan(sga.eps) & ~isnan(sga.pvl));
	for i=1:length(r)
		fprintf(fid, '%s\t%s\t%s\t%s\t%f\t%f\t%f\t%e\t%s\n',...
				sga.Cannon.Orf{r(i)},...
				sga.Cannon.Common{r(i)},...
				sga.Cannon.Orf{c(i)},...
				sga.Cannon.Common{c(i)},...
				sga.dbl(r(i),c(i)),...
				sga.eps(r(i),c(i)),...
				sga.dbl_std(r(i),c(i)),...
				sga.pvl(r(i),c(i)),...
				sga.sources{sga.src(r(i),c(i))});
	end
	fclose(fid);
end
