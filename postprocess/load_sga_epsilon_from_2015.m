function[sga] = load_sga_epsilon_from_2015(score_file_string, orf_file)
%function[sga] = load_sga_epsilon_from_release(score_file_string, orf_file)
   tic

   sga = struct()
   % make a Cannon
   Cannon = struct();
   Cannon.Orf = Csv2Cell(orf_file);
   Cannon.GENES = numel(Cannon.Orf);
   Cannon.Common = OrfToCommon(Cannon.Orf);
   Cannon.isArray = logical(zeros(1,Cannon.GENES));
   Cannon.isQuery = logical(zeros(Cannon.GENES, 1));

   Cannon.Map = java.util.HashMap(Cannon.GENES*2);
   Cannon.Map = Hash(Cannon.Map, Cannon.Orfs);
   Cannon.Map = Hash(Cannon.Map, Cannon.Common);

   % allocate matriciesk
   sga.eps = NaN(Cannon.GENES);
   sga.dbl = NaN(Cannon.GENES);
   sga.dbl_std = NaN(Cannon.GENES);
   sga.pvl = NaN(Cannon.GENES);
   sga.source= NaN(Cannon.GENES);

   sga.source_labels = {'DMA26', 'DMA30', 'TSA26', 'TSA30'}; % match previous where possible
   %valid_srcs = {'FG26', 'FG30', 'TS26', 'TS30', 'FG_MERGE', 'TS_MERGE'};
   src_map = Hash([], sga.source_labels);

   % read in the data
   % all numbers are doubles
   format = '%s%s%s%s%s%f%f%f%f%f%f';

   %burn the header
   fid = fopen(score_file_string, 'r');
   header = fgetl(fid);

   while ~feof(fid)
      segarray = textscan(fid, format, 1000, 'ReturnOnError', false);
      numlines = length(segarray{1});
      for i=1:numlines
         ixQ = Cannon.Map.get(segarray{1}{i});
         ixA = Cannon.Map.get(segarray{3}{i});

         if(isempty(ixQ) || isempty(ixA))
            fprintf('name map error %s %s\n', segarray{1}{i}, segarray{3}{i});
            continue;
         end

         % logical vectors
         Cannon.isQuery(ixQ) = true;
         Cannon.isArray(ixA) = true;

         % experiment
         sga.source(ixQ, ixA) = src_map.get(segarray{5}{i});

         % Epsilon Score =
         sga.eps(ixQ, ixA) = segarray{6}(i);

         % Pvalue pval
         sga.pvl(ixQ, ixA) = segarray{7}(i);

         % smfs
         sga.qfit(ixQ,ixA) = segarray{8}(i);
         sga.afit(ixQ,ixA) = segarray{9}(i);

         % Double mutant fitness dm_act
         sga.dbl(ixQ, ixA) = segarray{10}(i);
         sga.dbl_std(ixQ, ixA) = segarray{11}(i);

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

   sga.Cannon = Cannon;

   toc
end
