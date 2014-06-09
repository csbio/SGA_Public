function[] = MyCluster(mat, row_label1, row_label2, col_label, jobname, TYPE)
%function[] = MyCluster(mat, row_label1, row_label2, col_label, jobname, TYPE)
%
%	Defaults:
%	-e 9 0- inner product on arrays
%	-g 9 0- inner product on genes
%	-m a pairwise average linkage
% TYPE == 'inner_product' (default), or 'pearson'
%

	% default to inner_product for backward compatabiltity 
	if ~exist('TYPE', 'var')
		TYPE = 'inner_product';
	end

	if strcmp(TYPE, 'inner_product')
		command = 'cluster -e 9 -g 9 -m a -f '; % inner product
	elseif strcmp(TYPE, 'pearson')
		command = 'cluster -e 2 -g 2 -m a -f '; % pearson 
	else
		% give the user a chance to reset it and continue
		fprintf('invalid clustering method request (entering degubber)\n');
		keyboard
	end

	mat(isnan(mat)) = 0;
	filename = strcat(jobname, '.pcl');
	print_pcl_file(mat, row_label1, row_label2, col_label, filename);

	
	command = [command filename];
	disp(command)
	system(command);

	fix_tree_files(jobname);

	system(['rm -f ' filename]); % cleanup the pcl

end


function[] = fix_tree_files(jobname)
% this scales inner products to look like pearsons...

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
