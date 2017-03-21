function sga_ab = convert_bjv_to_ab(sga_bjv)
%function sga_ab = convert_bjv_to_ab(sga_bjv)
% extracts components of my SGA data structure to be 
% compatible with anastasisa's code

   sga_ab = struct();                                                                                                               
   sga_ab.queries = sga_bjv.Cannon.Orf(sga_bjv.Cannon.isQuery);                                                                             
   sga_ab.arrays  = sga_bjv.Cannon.Orf(sga_bjv.Cannon.isArray);                                                                             
   sga_ab.scores_eps  = sga_bjv.eps(sga_bjv.Cannon.isQuery, sga_bjv.Cannon.isArray);       
   sga_ab.scores_pvalue  = sga_bjv.pvl(sga_bjv.Cannon.isQuery, sga_bjv.Cannon.isArray);       

   sga_ab.queries_gn = StrainToAllele(sga_ab.queries);
   sga_ab.arrays_gn = StrainToAllele(sga_ab.arrays);

end