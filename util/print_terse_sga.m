function[] = print_terse_sga(sga, filename)
%function[] = print_terse_sga(sga, filename)
% Prints all non-nan scores from an sga struct
% Respects isQuery and isArray
% Orf Orf EPS PVL
	fid = fopen(filename, 'w');
	EPS = sga.eps(sga.Cannon.isQuery, sga.Cannon.isArray);
	PVL = sga.pvl(sga.Cannon.isQuery, sga.Cannon.isArray);

	% fprintf('using STRAIN names\n');
	% Qrf = sga.Cannon.Orf(sga.Cannon.isQuery);
	% Arf = sga.Cannon.Orf(sga.Cannon.isArray);

	fprintf('using ORF names\n');
	Qrf = StripOrfs(sga.Cannon.Orf(sga.Cannon.isQuery), 'first');
	Arf = StripOrfs(sga.Cannon.Orf(sga.Cannon.isArray), 'first');



	[r, c] = find(~isnan(EPS) & ~isnan(PVL));
	for i=1:length(r)
		fprintf(fid, '%s\t%s\t%f\t%e\n', Qrf{r(i)}, Arf{c(i)},EPS(r(i),c(i)),PVL(r(i),c(i)));
	end
end
