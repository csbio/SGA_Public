function[] = MyCluster(mat, row_label1, row_label2, col_label, jobname)
% MyCluster(mat, row_label1, row_label2, col_label, jobname)
%
%	Defaults:
%	-e 9 0- inner product on arrays
%	-g 9 0- inner product on genes
%	-m a pairwise average linkage
%

	mat(isnan(mat)) = 0;

	filename = strcat(jobname, '.pcl');
	print_pcl_file(mat, row_label1, row_label2, col_label, filename);

	command = 'cluster -e 9 -g 9 -m a -f ';
	command = [command filename];
	disp(command)
	system(command);


	fix_tree_files(jobname);

end


function[] = fix_tree_files(jobname)

%-% Array Tree
	filename = [jobname '.atr'];
	fid = fopen(filename, 'r');
	data = textscan(fid, '%s%s%s%f');
	fclose(fid);

	% shift then scale range -> [0 .. 1]
	data{4} = data{4} - min(data{4});
	data{4} = data{4} / max(data{4});

	% rewrite
	fid = fopen(filename, 'w');
	for i=1:length(data{1})
		fprintf(fid, '%s\t%s\t%s\t%f\n', ...
			data{1}{i},...
			data{2}{i},...
			data{3}{i},...
			data{4}(i));
	end
	fclose(fid);

%-% Gene Tree
	filename = [jobname '.gtr'];
	fid = fopen(filename, 'r');
	data = textscan(fid, '%s%s%s%f');
	fclose(fid);

	% shift then scale range -> [0 .. 1]
	data{4} = data{4} - min(data{4});
	data{4} = data{4} / max(data{4});

	% rewrite
	fid = fopen(filename, 'w');
	for i=1:length(data{1})
		fprintf(fid, '%s\t%s\t%s\t%f\n', ...
			data{1}{i},...
			data{2}{i},...
			data{3}{i},...
			data{4}(i));
	end
	fclose(fid);
end
