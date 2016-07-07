function[sga] = load_sga_epsilon_from_costanzo(score_file_string, orf_file)
%function[sga] = load_sga_epsilon_from_costanzo(score_file_string, orf_file)
tic

% make a Cannon
Cannon = struct();

fid = fopen(orf_file, 'r');
A = textscan(fid, '%s');
Cannon.Orf = A{1};
Cannon.Common = CommonToOrf(Cannon.Orf)
Cannon.GENES = length(Cannon.Orf);
Cannon.isArray = logical(zeros(1,Cannon.GENES));
Cannon.isQuery = logical(zeros(Cannon.GENES, 1));
fclose(fid);

Cannon.Map = hash_strings(Cannon.Orf);
Cannon.Map = hash_strings(Cannon.Common, Cannon.Map);

% Score file columns:
% Qorf Qcom Aorf Acom dm_act eps std pval experiment
% same as release
format = '%s%s%s%s%f%f%f%f%f%f%f%f%f'; % 4 strings, then all singles (64 bit pval)
block_size = 1000;
fid = fopen(score_file_string, 'r');

% allocate matriciesk
sga_eps = zeros(Cannon.GENES) + nan;
sga_dbl_std = zeros(Cannon.GENES) + nan;
sga_pvl = zeros(Cannon.GENES) + nan;


% Read block loop
while ~feof(fid)
	segarray = textscan(fid, format, block_size, 'ReturnOnError', false);
	numlines = length(segarray{1});
	for i=1:numlines
		ixQ = Cannon.Map.get(segarray{1}{i});
		ixA = Cannon.Map.get(segarray{3}{i});

		if(isempty(ixQ) || isempty(ixA))
			fprintf('name map error %s %s\n', segarray{1}{i}, segarray{3}{i});
			continue;
		end

		% Epsilon Score = 
		sga_eps(ixQ, ixA) = segarray{5}(i);

		% Double mutant fitness dm_act
		%sga_dbl(ixQ, ixA) = segarray{?}(i);
		sga_dbl_std(ixQ, ixA) = segarray{6}(i);
	

		% Pvalue pval
		sga_pvl(ixQ, ixA) = segarray{7}(i);

		% experiment...

		% boolean vectors
		% add these only if we saw a real value
		if(~isnan(sga_eps(ixQ, ixA)) && ~isnan(sga_pvl(ixQ, ixA)))
			Cannon.isArray(ixA) = true;
			Cannon.isQuery(ixQ) = true;
		end
	end
end
fclose(fid);

arrays = Cannon.Orf(Cannon.isArray);
Cannon.ArrayMap = hash_strings(arrays);

querys = Cannon.Orf(Cannon.isQuery);
Cannon.QueryMap = hash_strings(querys);

sga = struct('eps', sga_eps, 'pvl', sga_pvl, ...
   'Cannon', Cannon, 'dbl_std', sga_dbl_std);

toc
