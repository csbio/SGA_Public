function[sga] = load_sga_epsilon_from_release(score_file_string, orf_file)
%function[sga] = load_sga_epsilon_from_release(score_file_string)
tic

% make a Cannon
Cannon = struct();
Cannon.Map = java.util.HashMap(6000);

fid = fopen(orf_file, 'r');
A = textscan(fid, '%s');
Cannon.Orf = A{1};
Cannon.GENES = length(Cannon.Orf);
Cannon.isArray = boolean(zeros(1,Cannon.GENES));
Cannon.isQuery = boolean(zeros(Cannon.GENES, 1));
fclose(fid);

for i=1:Cannon.GENES
	Cannon.Map.put(java.lang.String(Cannon.Orf{i}), java.lang.Integer(i));
end

Cannon.Common = OrfToCommon(Cannon.Orf);
Cannon.Map = Hash(Cannon.Map, Cannon.Common);

% Score file columns:
% Qorf Qcom Aorf Acom dm_act eps std pval experiment
format = '%s%s%s%s%f32%f32%f32%f64%s'; % 4 strings, then all singles (64 bit pval)
block_size = 1000;
fid = fopen(score_file_string, 'r');
TOTAL_LINES = 14171105;
pct = 0;
iter_cnt = 0;


% allocate matriciesk
sga_eps = zeros(Cannon.GENES) + nan;
sga_dbl = zeros(Cannon.GENES) + nan;
sga_dbl_std = zeros(Cannon.GENES) + nan;
sga_pvl = zeros(Cannon.GENES) + nan;


% Read block loop
while ~feof(fid)
	segarray = textscan(fid, format, block_size, 'HeaderLines', 1, 'ReturnOnError', false);
	numlines = length(segarray{1});
	for i=1:numlines
		ixQ = Cannon.Map.get(segarray{1}{i});
		ixA = Cannon.Map.get(segarray{3}{i});

		if(isempty(ixQ) || isempty(ixA))
			fprintf('name map error %s %s\n', segarray{1}{i}, segarray{3}{i});
			continue;
		end

		% Epsilon Score = 
		sga_eps(ixQ, ixA) = segarray{6}(i);

		% Double mutant fitness dm_act
		sga_dbl(ixQ, ixA) = segarray{5}(i);
		sga_dbl_std(ixQ, ixA) = segarray{7}(i);
	

		% Pvalue pval
		sga_pvl(ixQ, ixA) = segarray{8}(i);

		% experiment...

		% boolean vectors
		Cannon.isArray(ixA) = true;
		Cannon.isQuery(ixQ) = true;


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
   'Cannon', Cannon, 'dbl_std', sga_dbl_std);

toc
