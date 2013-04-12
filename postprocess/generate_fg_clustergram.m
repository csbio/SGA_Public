function[] = generate_fg_clustergram(sga, name)
%function[] = generate_fg_clustergram(sga, name)


	fid = fopen('/project/csbio/lab_share/SGA/refdata/StrainID-Allele_map.csv', 'r');
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


	%{
	% CLUSTER GRAM
	% TECHS and ARRAY POSI
	a_labels = Orf_to_array_position(sga.Cannon.Orf(sga.Cannon.isArray),...
						sga.Cannon.allele(sga.Cannon.isArray));
	%q_labels = Add_Techs(sga.Cannon.Orf(sga.Cannon.isQuery),...
						%sga.Cannon.allele(sga.Cannon.isQuery));

	MyCluster(EPS,...
				 sga.Cannon.Orf(sga.Cannon.isQuery), ...
				 q_labels, ...
				 a_labels,...
				 [name '_CLUST']);
	%}

	% Zero insignif - data
	EPS(isnan(EPS)) = 0;
	EPS(EPS > -0.08 & EPS < 0.08) = 0;
	EPS(PVL > 0.05) = 0;
	
	MyCluster(EPS,...
				 sga.Cannon.Orf(sga.Cannon.isQuery), ...
				 sga.Cannon.allele(sga.Cannon.isQuery), ...
				 sga.Cannon.allele(sga.Cannon.isArray), [name '_CLUST']);
				 %sga.Cannon.Orf(sga.Cannon.isArray), [name '_CLUST']);
			 

	delete([name '_CLUST.pcl']);

	qorf = sga.Cannon.Orf(sga.Cannon.isQuery);
	qallele = allele(sga.Cannon.isQuery);

	aorf = sga.Cannon.Orf(sga.Cannon.isArray);
	aallele = allele(sga.Cannon.isArray);

	qix = ChromSort(qorf);
	aix = ChromSort(aorf);

	print_pcl_file(EPS(qix, aix), qallele(qix), qorf(qix), aorf(aix), [name '_CHROM.pcl']);
	%print_pcl_file(EPS(qix, aix), qallele(qix), qorf(qix), aallele(aix), [name '_CHROM.pcl']);


end


% Functions to add varying labels to a clustergram
%{
function[labels] = Orf_to_array_position(Orfs, Commons)

	fid = fopen('/project/csbio/lab_share/SGA/Main/postprocess/TS_Array_v6.csv', 'r');
	A = textscan(fid, '%s%s%s%s%s%s');
	A = [A(:)];
	fclose(fid);

	labels = cell(size(Orfs));
	for i=1:length(Orfs)
		ix = strmatch(upper(SStrain(Orfs{i})), A{3}, 'exact');
		if(isempty(ix))
			labels{i} = [Commons{i} ' ---'];
		else
			if(length(ix) > 4)
				Commons{i} = [Commons{i} ' ++ '	];
			end
			ix = ix(1); % grab the first one
			labels{i} = [Commons{i} sprintf(' P%s (%s %s)', A{4}{ix}, A{5}{ix}, A{6}{ix})];
		end
	end
end


function[sid] = SStrain(Orf)
	ix = strfind(Orf, '_');
	sid = Orf(ix+1:end);
end


function[labels] = Add_Techs(Orfs, Commons);

	fid = fopen('/project/csbio/lab_share/SGA/rawdata/120615/Techs.txt', 'r');
	A = textscan(fid, '%s%s%s'); % name date orf
	fclose(fid);

	labels = cell(size(Orfs));
	for i=1:length(Orfs)
		ix = strmatch(Orfs{i}, A{3});
		if(isempty(ix))
			labels{i} = [Commons{i} ' ??'];
		else
			ix = ix(1);
			labels{i} = [Commons{i} ' - ' A{1}{ix}];
		end
	end
end
%}
