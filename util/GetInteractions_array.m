function [orfs, coms, eps, pvl] = GetInteractions_array(sga, array, int_type, THRESH)
%function [orfs, coms, eps, pvl] = GetInteractions(sga, array, int_type='neg', THRESH=0.08)
% assumes the array exists, and is named exactly
% THRESH should be positive (abs())
% type defaults to 'neg' and THRESH to '0.08'

   if ~exist('int_type', 'var')
      int_type = 'neg';
   end

   if ~exist('THRESH', 'var')
      THRESH = 0.08;
   end

   assert(ismember(int_type, {'neg', 'pos'})); % TODO pos, both


   ixa = sga.Cannon.Map.get(array);
   if(isempty(ixa) || ~strcmp(sga.Cannon.Orf{ixa}, array))
      fprintf('warning Cannon.Map inconsistency... doing slow lookup\n')
      ixa = unique([strmatch(array, sga.Cannon.Common); ...
         strmatch(array, sga.Cannon.Common)])
   end

   if strcmp(int_type, 'neg')
      ixq = sga.eps(sga.Cannon.isQuery, ixa) <= -THRESH & sga.pvl(sga.Cannon.isQuery, ixa) < 0.05;
   elseif strcmp(int_type, 'pos')
      ixq = sga.eps(sga.Cannon.isQuery, ixa) >= THRESH & sga.pvl(sga.Cannon.isQuery, ixa) < 0.05;
   else
      fprintf('not implemented yet\n');
      return
   end

   queries = find(sga.Cannon.isQuery);
   orfs = sga.Cannon.Orf(queries(ixq));
   coms = sga.Cannon.Common(queries(ixq));
   eps = sga.eps(queries(ixq),ixa);
   pvl = sga.pvl(queries(ixq),ixa);
end