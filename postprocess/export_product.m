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
	



	dirname = split_by_delimiter('/', sga_outputfile);
	basename= split_by_delimiter('_', dirname{end});
	int_dirname = [join_by_delimiter(dirname(1:end-1), '/') '/interactions/'];
	system(['mkdir -p ' int_dirname]);


	% ------------------------ processed matfiles
		fields = split_by_delimiter('_', basename);
		project = fields{2};
		array   = fields{3};
		temp    = fields{4};
		construct = join_by_delimiter({project, array, temp}, '_');

	sga_raw = load_sga_epsilon_from_scorefile([sga_outputfile '.txt'], [sga_outputfile '.orf']);
	eval(sprintf('%s_raw = sga_raw;', construct)); % save a pre-filter version as well
	eval(sprintf('save %s%s_raw.mat %s_raw -v7.3', int_dirname, construct, construct)); % save to mat

	% layout file expects format (plate row col orf) [384] 16x24
	sga = filter_green_blocks_around_linkage(sga_raw, linkagefile, coord_file, layout_file, wild_type, border_strain);

	% this will do cobatch filter and AB BA disagreement filtereing in either case
	% but will not do any fitness depedant filtering unless 'false'
	[sga, fitness_struct]  = filter_interactions(sga, smfitnessfile, sga_inputfile, equiv_file);
   return

	% print out the cobatch filter receipt 
	cell2csv([sga_outputfile '_cobatch_removed.txt'], sga.Cannon.Orf(sga.cobatch_target));



	% ------------------------ clustergrams 
	clus_dirname = [join_by_delimiter(dirname(1:end-1), '/') '/clustergrams/'];
	system(['mkdir -p ' clus_dirname]);
	clus_basename = basename;
	clus_basename{1} = 'clustergram';
	clus_basename = join_by_delimiter(clus_basename, '_');

   % skip the clustergrams
	generate_fg_clustergram(sga, [clus_dirname clus_basename]);


	% ------------------------ interactions / profiles

	% int_basename{1} = 'release';
	% prof_basename = basename;
	% prof_basename{1} = 'profile';
	% prof_basename = join_by_delimiter(prof_basename, '_');
	% print_profile_data(sga, [int_dirname '/' prof_basename '.txt']);


	eval(sprintf('%s = sga;', construct)); % rename struct
	eval(sprintf('save %s%s.mat %s  -v7.3', int_dirname, construct, construct)); % save to mat


end



function[] = print_profile_data(sga, outputfile);
	%{
	----------------------------------------------------------
	TYPE X?
	OUTPUT; RELEASE SGAdata
	9 col
	----------------------------------------------------------
	1  Query ORF 
	2  Query gene name
	3  Array ORF 
	4  Array gene name
	5  Genetic interaction score (eps)
	6  Standard deviation
	7  p-value
	8  Double mutant fitness
	9  boolean intermediate cutoff
	%}

	fid = fopen(outputfile, 'w');


	[r, c] = find(~isnan(sga.eps) & ~isnan(sga.pvl));
	for i=1:length(r)
		if(abs(sga.eps(r(i),c(i))) > 0.08 && sga.pvl(r(i),c(i)) < 0.05)
			int_status = 'TRUE';
		else
			int_status = 'FALSE';
		end

		fprintf(fid, '%s\t%s\t%s\t%s\t%f\t%f\t%e\t%f\t%s\n',...
				sga.Cannon.Orf{r(i)},...
				sga.Cannon.Common{r(i)},...
				sga.Cannon.Orf{c(i)},...
				sga.Cannon.Common{c(i)},...
				sga.eps(r(i),c(i)),...
				sga.dbl_std(r(i),c(i)),...
				sga.pvl(r(i),c(i)),...
				sga.dbl(r(i),c(i)),...
				int_status);
	end
	fclose(fid);
end
