function[sga] = load_sga_epsilon_from_2016_release(score_file_string)
%function[sga] = load_sga_epsilon_from_2016_release(score_file_string)
tic

strain_file = sprintf('%s_xxx_strainlist.tmp', score_file_string);
% generate a strainlist file in the shell
exec_str1 = sprintf('system(''cut -f1 %s | sort -u >  %s'')', score_file_string, strain_file);
exec_str2 = sprintf('system(''cut -f3 %s | sort -u >> %s'')', score_file_string, strain_file);

retval = eval(exec_str1);
if retval == 0
   retval = eval(exec_str2);
end

if retval > 0
   error('trouble making strainlist file');
end


% make a Cannon
Cannon = struct();
Cannon.Map = java.util.HashMap(10000);

fid = fopen(strain_file, 'r');
A = textscan(fid, '%s');
fclose(fid);
system(sprintf('rm %s', strain_file));

Cannon.Orf = A{1};
Cannon.Common = cell(size(Cannon.Orf));
Cannon.GENES = length(Cannon.Orf);

for i=1:Cannon.GENES
	Cannon.Map.put(java.lang.String(Cannon.Orf{i}), java.lang.Integer(i));
end

Cannon.sources = {'DMA26', 'DMA30', 'TSA26', 'TSA30'};
Cannon.src_map = Hash([], Cannon.sources);

% Score file columns:
% Qstrain Qallele Astrain Aallele Arraytype/Temp epsilon pval qsmf asmf dmf dmf_std
format = '%s%s%s%s%s%f%f%f%f%f%f'; 
block_size = 1000;
fid = fopen(score_file_string, 'r');

% allocate matriciesk
sga_eps = zeros(Cannon.GENES) + nan;
sga_pvl = zeros(Cannon.GENES) + nan;
sga_dbl = zeros(Cannon.GENES) + nan;
sga_dbl_std = zeros(Cannon.GENES) + nan;
sga_src = zeros(Cannon.GENES) + nan;
sga_qfit = zeros(Cannon.GENES) + nan;
sga_afit = zeros(Cannon.GENES) + nan;

% toss the header
header = fgetl(fid);
if(~strcmp(header(1:5), 'Query')) % no header, rewind
   frewind(fid);
end

% Read block loop
while ~feof(fid)
	segarray = textscan(fid, format, block_size, 'ReturnOnError', false, 'Delimiter', '\t');
	numlines = length(segarray{1});
	for i=1:numlines
		ixQ = Cannon.Map.get(segarray{1}{i});
		ixA = Cannon.Map.get(segarray{3}{i});

      % fill in common (allele) if we haven't seen it yet
      if(isempty(Cannon.Common{ixQ}))
         Cannon.Common{ixQ} = segarray{2}{i};
      end
      if(isempty(Cannon.Common{ixA}))
         Cannon.Common{ixA} = segarray{4}{i};
      end

      % Source info
      sga_src(ixQ, ixA) = Cannon.src_map.get(segarray{5}{i});

		% Epsilon Score  
		sga_eps(ixQ, ixA) = segarray{6}(i);

		% Pvalue pval
		sga_pvl(ixQ, ixA) = segarray{7}(i);

      % SMF
      sga_qfit(ixQ, ixA) = segarray{8}(i);
      sga_afit(ixQ, ixA) = segarray{9}(i);

		% Double mutant fitness
		sga_dbl(ixQ, ixA) = segarray{10}(i);
		sga_dbl_std(ixQ, ixA) = segarray{11}(i);
	

	end
end
fclose(fid);

% boolean vectors
Cannon.isQuery = sum(~isnan(sga_eps), 2) > 0;
Cannon.isArray = sum(~isnan(sga_eps), 1) > 0;

sga = struct('eps', sga_eps, 'pvl', sga_pvl, ...
   'Cannon', Cannon, 'dbl', sga_dbl, 'dbl_std', sga_dbl_std, ...
   'src', sga_src, 'qfit', sga_qfit, 'afit', sga_afit);

toc
