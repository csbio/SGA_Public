function[sga] = load_sga_epsilon_from_scorefile(score_file_string, unique_names_file)
%function[sga] = load_sga_epsilon_from_scorefile(score_file_string, unique_names_file)
tic

% make a Cannon
Cannon = struct();
if(~exist('unique_names_file', 'var'))
	unique_names_file = '/project/csbio/lab_share/SGA/Triple/scored/genelist_unique_Y.txt';
end

fid = fopen(unique_names_file, 'r');
A = textscan(fid, '%s');
Cannon.Orf = A{1};
Cannon.GENES = length(Cannon.Orf);
Cannon.isArray = false(1,Cannon.GENES);
Cannon.isQuery = false(Cannon.GENES, 1);
fclose(fid);

Cannon.Common = OrfToCommon(Cannon.Orf);
Cannon.Map = hash_strings(Cannon.Orf);
Cannon.Map = hash_strings(Cannon.Common, Cannon.Map);

% Score file columns:
% Qorf Aorf escore std pval smfit1 std smfit2 std dm_exp dm_act std
format = '%s%s%f32%f32%f64%f32%f32%f32%f32%f32%f32%f32'; % 2 strings, then all singles
block_size = 1000;
fid = fopen(score_file_string, 'r');

% allocate matriciesk
sga_eps = zeros(Cannon.GENES) + nan;
sga_dbl = zeros(Cannon.GENES) + nan;
sga_dbl_std = zeros(Cannon.GENES) + nan;
sga_escore = zeros(Cannon.GENES) + nan;
sga_pvl = zeros(Cannon.GENES) + nan;
sga_fit = zeros(Cannon.GENES,1) +nan;


% Read block loop
while ~feof(fid)
	segarray = textscan(fid, format, block_size, 'ReturnOnError', false);
	numlines = length(segarray{1});
	for i=1:numlines
		ixQ = Cannon.Map.get(segarray{1}{i});
		ixA = Cannon.Map.get(segarray{2}{i});

		if(isempty(ixQ) || isempty(ixA))
			fprintf('name map error %s %s\n', segarray{1}{i}, segarray{2}{i});
			continue;
		end

		% Epsilon Score = dm_act - dm_exp
		sga_eps(ixQ, ixA) = segarray{11}(i) - segarray{10}(i);

		% Double mutant fitness dm_act
		sga_dbl(ixQ, ixA) = segarray{11}(i);
		sga_dbl_std(ixQ, ixA) = segarray{12}(i);
	

		% Pvalue pval
		sga_pvl(ixQ, ixA) = segarray{5}(i);

		% boolean vectors
		Cannon.isArray(ixA) = true;
		Cannon.isQuery(ixQ) = true;

		% fitness
		sga_fit(ixQ) = segarray{6}(i);
		sga_fit(ixA) = segarray{8}(i);
		
		sga_escore(ixQ, ixA) = segarray{3}(i);
		

	end
end
fclose(fid);

arrays = Cannon.Orf(Cannon.isArray);
Cannon.ArrayMap = hash_strings(arrays);

querys = Cannon.Orf(Cannon.isQuery);
Cannon.QueryMap = hash_strings(querys);

sga = struct('eps', sga_eps, 'dbl', sga_dbl, 'pvl', sga_pvl, ...
   'Cannon', Cannon, 'dbl_std', sga_dbl_std, 'escore', sga_escore, ...
	'fit', sga_fit);

toc
