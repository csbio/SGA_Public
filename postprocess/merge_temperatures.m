function[sga] = merge_temperatures(primary, p_name, secondary, s_name)
%function[sga] = merge_temperatures(primary, p_name, secondary, s_name)
% takes in two structs and labels, and produces a merged struct.
% respects isQuery/isArray
% values from PRIMARY override SECONDARY
% all valid interactions will survive (union)
% labels for each interaction source will get saved

   % make a new cannon
   % use the intersection of arrays
   arrays = intersect(primary.Cannon.Orf(primary.Cannon.isArray),...
                  secondary.Cannon.Orf(secondary.Cannon.isArray));

   % unset arrays not in the intersection
   primary.Cannon.isArray(~ismember(primary.Cannon.Orf, arrays)) = false;
   secondary.Cannon.isArray(~ismember(secondary.Cannon.Orf, arrays)) = false;

   queries = union(primary.Cannon.Orf(primary.Cannon.isQuery),...
                  secondary.Cannon.Orf(secondary.Cannon.isQuery));

   sga = struct();
   sga.Cannon = struct();

   sga.Cannon.Orf = [queries; arrays];
   sga.Cannon.Common = OrfToCommon(sga.Cannon.Orf);
   sga.Cannon.GENES = length(sga.Cannon.Orf);

   sga.Cannon.isQuery = false(sga.Cannon.GENES,1);
   sga.Cannon.isQuery(1:length(queries)) = true;
   sga.Cannon.isArray = false(1,sga.Cannon.GENES);
   sga.Cannon.isArray(end-length(arrays)+1:end) = true;

   sga.Cannon.Map = Hash([], sga.Cannon.Orf);
   sga.Cannon.QueryMap = Hash([], queries);
   sga.Cannon.ArrayMap = Hash([], arrays);
   
   sga.eps = NaN(sga.Cannon.GENES);
   sga.dbl = NaN(sga.Cannon.GENES);
   sga.pvl = NaN(sga.Cannon.GENES);
   sga.dbl_std = NaN(sga.Cannon.GENES);
   sga.escore = NaN(sga.Cannon.GENES);
   sga.qfit = NaN(sga.Cannon.GENES);
   sga.afit = NaN(sga.Cannon.GENES);
   sga.source = NaN(sga.Cannon.GENES);
   sga.source_labels = {p_name, s_name};

   % put all valid secondary values in first
   % tmp variables in pre-merge orientation
   EPS = secondary.eps(secondary.Cannon.isQuery, secondary.Cannon.isArray);
   DBL = secondary.dbl(secondary.Cannon.isQuery, secondary.Cannon.isArray);
   PVL = secondary.pvl(secondary.Cannon.isQuery, secondary.Cannon.isArray);
   DBL_STD = secondary.dbl_std(secondary.Cannon.isQuery, secondary.Cannon.isArray);
   ESCORE = secondary.escore(secondary.Cannon.isQuery, secondary.Cannon.isArray);
   QFIT = repmat(secondary.fit(secondary.Cannon.isQuery), 1, sum(secondary.Cannon.isArray));
   AFIT = repmat(secondary.fit(secondary.Cannon.isArray)', sum(secondary.Cannon.isQuery), 1);

   % map destination in absolute coords
   qix = apply_map(sga.Cannon.Map, secondary.Cannon.Orf(secondary.Cannon.isQuery));
   aix = apply_map(sga.Cannon.Map, secondary.Cannon.Orf(secondary.Cannon.isArray));

   sga.eps(qix, aix) = EPS;
   sga.dbl(qix, aix) = DBL;
   sga.pvl(qix, aix) = PVL;
   sga.dbl_std(qix, aix) = DBL_STD;
   sga.escore(qix, aix) = ESCORE;
   sga.qfit(qix,aix) = QFIT;
   sga.afit(qix, aix) = AFIT;
   sga.source(qix, aix) = ones(size(EPS))+1;


   % primary overwrite -------
   % since array sets are identical, we take all values
   % in effect, replacing entire rows

   EPS = primary.eps(primary.Cannon.isQuery, primary.Cannon.isArray);
   DBL = primary.dbl(primary.Cannon.isQuery, primary.Cannon.isArray);
   PVL = primary.pvl(primary.Cannon.isQuery, primary.Cannon.isArray);
   DBL_STD = primary.dbl_std(primary.Cannon.isQuery, primary.Cannon.isArray);
   ESCORE = primary.escore(primary.Cannon.isQuery, primary.Cannon.isArray);
   QFIT = repmat(primary.fit(primary.Cannon.isQuery), 1, sum(primary.Cannon.isArray));
   AFIT = repmat(primary.fit(primary.Cannon.isArray)', sum(primary.Cannon.isQuery), 1);

   % map destination in absolute coords
   qix = apply_map(sga.Cannon.Map, primary.Cannon.Orf(primary.Cannon.isQuery));
   aix = apply_map(sga.Cannon.Map, primary.Cannon.Orf(primary.Cannon.isArray));

   sga.eps(qix, aix) = EPS;
   sga.dbl(qix, aix) = DBL;
   sga.pvl(qix, aix) = PVL;
   sga.dbl_std(qix, aix) = DBL_STD;
   sga.escore(qix, aix) = ESCORE;
   sga.qfit(qix, aix) = QFIT;
   sga.afit(qix, aix) = AFIT;
   sga.source(qix,aix) = 1;

end
%{
trip_fg_t30_raw = 

        eps: [8076x8076 double]
        dbl: [8076x8076 double]
        pvl: [8076x8076 double]
     Cannon: [1x1 struct]
    dbl_std: [8076x8076 double]
     escore: [8076x8076 double]
        fit: [8076x1 double]

>> trip_fg_t30_raw.Cannon

ans = 

         Map: [12036 java.util.HashMap]
         Orf: {8076x1 cell}
       GENES: 8076
     isArray: [1x8076 logical]
     isQuery: [8076x1 logical]
      Common: {8076x1 cell}
    ArrayMap: [3841 java.util.HashMap]
    QueryMap: [3780 java.util.HashMap]
%}
function [ids] = apply_map(map, strings)
   ids = zeros(size(strings));
   for i=1:numel(strings)
      ix = map.get(strings{i});
      if ~isempty(ix)
         ids(i) = ix;
      end
   end
end
   
