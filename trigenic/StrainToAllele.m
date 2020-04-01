function[result, allele_map] = StrainToAllele(cellarr, mappingfile)
%function[result] = StrainToAllele(cellarr, [mappingfile | allele_map])
   
   if(~exist('mappingfile', 'var') || isempty(mappingfile))
      mappingfile = '~/Research/Data/YeastGeneMap/strain_orf_common_allele_160425.txt';
      fprintf('using default mapping file: %s\n', mappingfile);
   end

  if ischar(mappingfile)
     fid = fopen(mappingfile, 'r');
     A = textscan(fid, '%s%s%s%s', 'Delimiter', '\t', 'ReturnOnError', false);
     fclose(fid);
     allele_map = [A{1} A{4}];
   elseif iscell(mappingfile)
      allele_map = mappingfile;
   end
   map = Hash([], allele_map(:,1));

   result = cellarr; % seed with default values
   for i=1:numel(cellarr)
      ix = map.get(cellarr{i});
      if(~isempty(ix))
         result{i} = allele_map{ix,2};

      % if not found, and is dma or sn, strip the tag and lower the orf
      elseif(length(findstr('_dma', cellarr{i})) > 0)
         result{i} = lower(cellarr{i}(1:findstr('_dma', cellarr{i})-1));

      elseif(length(findstr('_sn', cellarr{i})) > 0)
         result{i} = lower(cellarr{i}(1:findstr('_sn', cellarr{i})-1));

      elseif(length(findstr('_nan', cellarr{i})) > 0)
         result{i} = lower(cellarr{i}(1:findstr('_nan', cellarr{i})-1));

      end
   end
end

