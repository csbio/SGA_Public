function[] = generate_fg_clustergram(sga, name)
%function[] = generate_fg_clustergram(sga, name)
% produces a CLUST and CHROM version


	%fid = fopen('/project/csbio/lab_share/SGA/refdata/StrainID-Allele_map_140107.csv', 'r');
	fid = fopen('/project/csbio/benjamin/Data/StrainID-Allele_map_140116.csv', 'r');

	A = textscan(fid, '%s%s');
	fclose(fid);

	m = java.util.HashMap();
	for i=1:1:length(A{1})
		m.put(java.lang.String(A{1}{i}), java.lang.String(A{2}{i}));
	end



	allele = sga.Cannon.Orf;
	for i=1:length(allele)
		ix = strfind(allele{i}, '_ts');
		if isempty(ix)
			allele{i} = sga.Cannon.Common{i};
			continue
		end

		a = m.get(allele{i}(ix+1:end));
		if isempty(a)
			continue
		end

		allele{i} = a;
	end
	sga.Cannon.allele = allele;

	EPS = sga.eps(sga.Cannon.isQuery, sga.Cannon.isArray);
	PVL = sga.pvl(sga.Cannon.isQuery, sga.Cannon.isArray);

	qorf = sga.Cannon.Orf(sga.Cannon.isQuery);
	qallele = allele(sga.Cannon.isQuery);
	aorf = sga.Cannon.Orf(sga.Cannon.isArray);
	aallele = allele(sga.Cannon.isArray);

	qix = ChromSort(qorf);
	aix = ChromSort(aorf);

	print_pcl_file(EPS(qix, aix), qallele(qix), qorf(qix), aorf(aix), [name '_CHROM.pcl']);
	%print_pcl_file(EPS(qix, aix), qallele(qix), qorf(qix), aallele(aix), [name '_CHROM.pcl']);


	% Zero insignif - data
	fprintf('Zeroing insignificant data\n');
	EPS(isnan(EPS)) = 0;
	EPS(EPS > -0.08 & EPS < 0.08) = 0;
	EPS(PVL > 0.05) = 0;
	%fprintf('NOT  Zeroing insignificant data\n');
	
	MyCluster(EPS,...
				 sga.Cannon.Orf(sga.Cannon.isQuery), ...
				 sga.Cannon.allele(sga.Cannon.isQuery), ...
				 sga.Cannon.allele(sga.Cannon.isArray), [name '_CLUST']);
				 %sga.Cannon.Orf(sga.Cannon.isArray), [name '_CLUST']);
			 


end
