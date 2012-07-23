function[sga] = export_product(sga_inputfile, sga_outputfile, smfitnessfile, ...
              linkagefile, coord_file, layout_file, equiv_file, wild_type, border_strain)
%function[] = export_product(sga_inputfile, sga_outputfile, smfitnessfile, ...
%                  linkagefile, coord_file, layout_file, wild_type, border_strain)
% post processing pipeline...
% Loads the datafile (it is stupid to write the file, then load it, but at least it keeps things modular)
% green-block
% interaction filters the data
% prints interaction data
% prints profile data
% makes clustergrams

% sga_output file = outputfile from compute_sgascore
% ie, it is extensionless
	
	sga_raw = load_sga_epsilon_from_scorefile([sga_outputfile '.txt'], [sga_outputfile '.orf']);

	% layout file expects format (plate row col orf) [384] 16x24
	sga = filter_green_blocks_around_linkage(sga_raw, linkagefile, coord_file, layout_file, wild_type, border_strain);

	[sga, fitness_struct]  = filter_interactions(sga, smfitnessfile, sga_inputfile, equiv_file);

	% ------------------------ clustergrams 
	dirname = split_by_delimiter('/', sga_outputfile);
	basename= split_by_delimiter('_', dirname{end});
	
	clus_dirname = [join_by_delimiter(dirname(1:end-1), '/') '/clustergrams/'];
	clus_basename = basename;
	clus_basename{1} = 'clustergram';
	clus_basename = join_by_delimiter(clus_basename, '_');
	system(['mkdir -p ' clus_dirname]);
	generate_fg_clustergram(sga, [clus_dirname clus_basename]);




	% ------------------------ interactions / profiles

	int_dirname = [join_by_delimiter(dirname(1:end-1), '/') '/interactions/'];
	int_basename = basename;
	int_basename{1} = 'interaction';
	system(['mkdir -p ' int_dirname]);
	int_basename = join_by_delimiter(int_basename, '_');
	print_interaction_data(sga, fitness_struct, [int_dirname '/' int_basename '.txt']);


	prof_basename = basename;
	prof_basename{1} = 'profile';
	prof_basename = join_by_delimiter(prof_basename, '_');
	print_profile_data(sga, fitness_struct, [int_dirname '/' prof_basename '.txt']);

	


	% ------------------------ processed matfiles
		fields = split_by_delimiter('_', basename);
		project = fields{2};
		array   = fields{3};
		temp    = fields{4};
		construct = join_by_delimiter({project, array, temp}, '_');

	eval(sprintf('%s = sga;', construct)); % rename struct
	eval(sprintf('save %s%s.mat %s', int_dirname, construct, construct)); % save to mat

	eval(sprintf('%s_raw = sga_raw;', construct)); % save a pre-filter version as well
	eval(sprintf('save %s%s_raw.mat %s_raw', int_dirname, construct, construct)); % save to mat

end


function[] = print_interaction_data(sga, fitness_struct, outputfile);
	%{
	for now we'll match the old format...
	------
	TYPE 2
	OUTPUT; published SGAdata
	13 col: example sgadata_costanzo2009_rawdata_101120.txt.gz
	----------------------------------------------------------
	1  Query ORF 
	2  Query gene name
	3  Array ORF 
	4  Array gene name
	5  Genetic interaction score (eps)
	6  Standard deviation
	7  p-value
	8  Query single mutant fitness (SMF)
	9  Query SMF standard deviation
	10 Array SMF 
	11 Array SMF standard deviation
	12 Double mutant fitness
	13 Double mutant fitness standard deviation 
	%}


	% re-arrange fitness data according to cannon
	fitness = nan(sga.Cannon.GENES, 2);
	for i=1:sga.Cannon.GENES
		ix = strmatch(sga.Cannon.Orf{i}, fitness_struct(:,1), 'exact');
		if(~isempty(ix))
			fitness(i,:) = cell2mat(fitness_struct(ix,[2,3]));
		end
	end

	fid = fopen(outputfile, 'w');

	epsilon = sga.eps;
	epsilon(isnan(epsilon)) = 0;
	epsilon(epsilon > -0.08 & epsilon < 0.08) = 0;

	[r, c] = find(epsilon);
	for i=1:length(r)
		fprintf(fid, '%s\t%s\t%s\t%s\t%f\t%f\t%e\t%f\t%f\t%f\t%f\t%f\t%f\n',...
				sga.Cannon.Orf{r(i)},...
				sga.Cannon.Common{r(i)},...
				sga.Cannon.Orf{c(i)},...
				sga.Cannon.Common{c(i)},...
				sga.eps(r(i),c(i)),...
				sga.dbl_std(r(i),c(i)),...
				sga.pvl(r(i),c(i)),...
				fitness(r(i),1),...
				fitness(r(i),2),...
				fitness(c(i),1),...
				fitness(c(i),2),...
				sga.dbl(r(i),c(i)),...
				sga.dbl_std(r(i),c(i)));
	end
	fclose(fid);
end

function[] = print_profile_data(sga, fitness_struct, outputfile);
	% same as above but insignificant interactions are printed also

	% re-arrange fitness data according to cannon
	fitness = nan(sga.Cannon.GENES, 2);
	for i=1:sga.Cannon.GENES
		ix = strmatch(sga.Cannon.Orf{i}, fitness_struct(:,1), 'exact');
		if(~isempty(ix))
			fitness(i,:) = cell2mat(fitness_struct(ix,[2,3]));
		end
	end

	fid = fopen(outputfile, 'w');

	epsilon = sga.eps;
	epsilon(isnan(epsilon)) = 0;
%	epsilon(epsilon > -0.08 & epsilon < 0.08) = 0;

	[r, c] = find(epsilon);
	for i=1:length(r)
		fprintf(fid, '%s\t%s\t%s\t%s\t%f\t%f\t%e\t%f\t%f\t%f\t%f\t%f\t%f\n',...
				sga.Cannon.Orf{r(i)},...
				sga.Cannon.Common{r(i)},...
				sga.Cannon.Orf{c(i)},...
				sga.Cannon.Common{c(i)},...
				sga.eps(r(i),c(i)),...
				sga.dbl_std(r(i),c(i)),...
				sga.pvl(r(i),c(i)),...
				fitness(r(i),1),...
				fitness(r(i),2),...
				fitness(c(i),1),...
				fitness(c(i),2),...
				sga.dbl(r(i),c(i)),...
				sga.dbl_std(r(i),c(i)));
	end
	fclose(fid);
end
