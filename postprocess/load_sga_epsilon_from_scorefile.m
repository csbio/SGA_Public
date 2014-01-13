function[sga] = load_sga_epsilon_from_scorefile(score_file_string, unique_names_file)
%function[sga] = load_sga_epsilon_from_scorefile(score_file_string, unique_names_file)
tic

% make a Cannon
Cannon = struct();
if(~exist('unique_names_file', 'var'))
	unique_names_file = '/project/csbio/lab_share/SGA/Triple/scored/genelist_unique_Y.txt';
end

fid = fopen(unique_names_file, 'r');
Cannon.Map = java.util.HashMap(6000);
A = textscan(fid, '%s');
Cannon.Orf = A{1};
Cannon.GENES = length(Cannon.Orf);
Cannon.isArray = logical(zeros(1,Cannon.GENES));
Cannon.isQuery = logical(zeros(Cannon.GENES, 1));
fclose(fid);

for i=1:Cannon.GENES
	%Cannon.Map.put(Cannon.Orf{i}, i);
	Cannon.Map.put(java.lang.String(Cannon.Orf{i}), java.lang.Integer(i));
end

Cannon.Common = OrfToCommon(Cannon.Orf);
Cannon.Map = Hash(Cannon.Map, Cannon.Common);

% Score file columns:
% Qorf Aorf escore std pval smfit1 std smfit2 std dm_exp dm_act std
format = '%s%s%f32%f32%f64%f32%f32%f32%f32%f32%f32%f32'; % 2 strings, then all singles
block_size = 1000;
fid = fopen(score_file_string, 'r');
TOTAL_LINES = 14171105;
pct = 0;
iter_cnt = 0;


% allocate matriciesk
sga_eps = zeros(Cannon.GENES) + nan;
sga_dbl = zeros(Cannon.GENES) + nan;
sga_dbl_std = zeros(Cannon.GENES) + nan;
sga_escore = zeros(Cannon.GENES) + nan;
sga_pvl = zeros(Cannon.GENES) + nan;
sga_fit = zeros(Cannon.GENES,1) +nan;


% Read block loop
while ~feof(fid)
	segarray = textscan(fid, format, block_size);
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
Cannon.ArrayMap = java.util.HashMap(length(arrays));
for i=1:length(arrays)
	Cannon.ArrayMap.put(arrays{i}, i);
end

querys = Cannon.Orf(Cannon.isQuery);
Cannon.QueryMap = java.util.HashMap(length(querys));
for i=1:length(querys)
	Cannon.QueryMap.put(querys{i}, i);
end

sga = struct('eps', sga_eps, 'dbl', sga_dbl, 'pvl', sga_pvl, ...
   'Cannon', Cannon, 'dbl_std', sga_dbl_std, 'escore', sga_escore, ...
	'fit', sga_fit);

toc
