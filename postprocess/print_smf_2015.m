function [table] = print_smf_2015(file_26, file_30, fg_merge, ts_merge)

   all_strains = fg_merge.Cannon.Orf(fg_merge.Cannon.isQuery | fg_merge.Cannon.isArray');
   all_strains = union(all_strains, ts_merge.Cannon.Orf(ts_merge.Cannon.isQuery | ts_merge.Cannon.isArray'));
   all_orfs = StripOrfs(all_strains);
   all_alleles = StrainToAllele(all_strains);

   % all the strains in the dataset will form the index
   table_header = {'StrainID', 'ORF', 'AlleleName', 'smf20°', 'smf20°-stddev', 'smf30°', 'smf30°-stddev'};

   table = cell(length(all_strains), length(table_header));
   table(:,1) = all_strains;
   table(:,2) = all_orfs;
   table(:,3) = all_alleles;

   % fill in 26 
   dat26 = Csv2Cell(file_26);
   h26 = hash_strings(dat26(:,1));
   for i=1:length(all_strains)
      ix = h26.get(all_strains{i});
      if isempty(ix)
         table(i,4:5) = {'NaN', 'NaN'};
      else
         table(i,4:5) = dat26(i,2:3);
      end
   end

   % fill in 30
   dat30 = Csv2Cell(file_30);
   h30 = hash_strings(dat30(:,1));
   for i=1:length(all_strains)
      ix = h30.get(all_strains{i});
      if isempty(ix)
         table(i,6:7) = {'NaN', 'NaN'};
      else
         table(i,6:7) = dat30(i,2:3);
      end
   end

   table = [table_header; table];
end

