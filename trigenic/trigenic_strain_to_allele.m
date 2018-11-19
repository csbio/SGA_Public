function [alleles] = trigenic_strain_to_allele(strains)
% function [alleles] = trigenic_strain_to_allele(strains)
%
% Calls StrainToAllele

   % works for single string
   if isstr(strains)
      strains = {strains};
   end


   % map alleles for normal (array?) strains
   [~, DM_strain_ix] = substrmatch('+', strains); % boolean

   % don't try to map / convert controls
   URA_ix = strcmp('URA3control+YDL227C_y13096', strains);
   DM_strain_ix(URA_ix) = false;

   alleles = cell(size(strains));
   alleles(~DM_strain_ix) = StrainToAllele(strains(~DM_strain_ix));

   % Build a strain table for TM data
   trigenic_strain_file = 'assignment_file_170328.csv';
   fid = fopen(trigenic_strain_file, 'r');
   %StrainID1       ORF1    StrainID2       ORF2    DMstrainID      SM1strainID     SM2strainID     Annotation
   A = textscan(fid, repmat('%s', 1, 8), 'ReturnOnError', false, 'Delimiter', '\t', 'HeaderLines', 1);
   fclose(fid);

   % enforce some conventions on case
   A{1} = lower(A{1});
   A{2} = upper(A{2});
   A{3} = lower(A{3});
   A{4} = upper(A{4});
   A{5} = lower(A{5});
   A{6} = lower(A{6});
   A{7} = lower(A{7});
   A{8} = upper(A{8});

   % convert to a straight 2D cell array
   TS = A{1};
   for i=2:8
      TS = [TS A{i}]; 
   end

   % build a table of all SM tm-ids to strain ids -> alleles
   tm_s1 = [JoinOrfs(TS(:,2), TS(:,6), '_'), JoinOrfs(TS(:,2), TS(:,1), '_')]; % [SM1_tm##   SM1_dma##] (e.g.)
   tm_s2 = [JoinOrfs(TS(:,4), TS(:,7), '_'), JoinOrfs(TS(:,4), TS(:,3), '_')]; % [SM2_tm##   SM2_dma##]
   
   % run a normal allele lookup on source strains
   tm_al1 = StrainToAllele(tm_s1(:,2)); 
   tm_al2 = StrainToAllele(tm_s2(:,2));

   dm_s = JoinOrfs(JoinOrfs(TS(:,2), TS(:,4), '+'), TS(:,5), '_'); % X+Y_tm##
   dm_al = JoinOrfs(tm_al1, tm_al2, '+');                          % allele(x)+allele(y)

   TM_strain_table = [[tm_s1(:,1) tm_al1]; [tm_s2(:,1) tm_al2]; [dm_s dm_al]];

   map = hash_strings(TM_strain_table(:,1));
   ix = apply_map(map, StripHO(strains(DM_strain_ix)));
   alleles(DM_strain_ix) = TM_strain_table(ix,2);


