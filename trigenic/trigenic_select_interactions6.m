function [ map_eps, trip ] = trigenic_select_interactions6(trip)
   % this is a fork of 5 that adds back support of positive interaction classes for paralog paper
   % returns a map in the global space with trigenic interaction classification
   % Classes 
   % [+ | -][Novel | [[ Q+ | Q- ] [ A+ | A- | A+- ]]]


   thresh = 0.08;

   % annotate all base interactions in place
   map_eps = cell(size(trip.eps));
   map_eps(:) = {'x'}; 

   neg_int_all = (trip.eps < -thresh) & (trip.pvl < 0.05);
   map_eps(neg_int_all) = {'-'};

   pos_int_all = (trip.eps >  thresh) & (trip.pvl < 0.05);
   map_eps(pos_int_all) = {'+'};

   % prepend D or T for digen and trigen
   [~, ~, ix_all] = find_trigenic_players(trip);
   map_eps(ix_all(:,1), trip.Cannon.isArray) = cellfun(@prepT, map_eps(ix_all(:,1), trip.Cannon.isArray), 'UniformOutput', false);
   map_eps(ix_all(:,2), trip.Cannon.isArray) = cellfun(@prepD, map_eps(ix_all(:,2), trip.Cannon.isArray), 'UniformOutput', false);
   map_eps(ix_all(:,3), trip.Cannon.isArray) = cellfun(@prepD, map_eps(ix_all(:,3), trip.Cannon.isArray), 'UniformOutput', false);
   % digen are now done, but we still need the interactions to classify trigenics

   % extract from struct
   % mode of dmq only (novel / modified) ======================================================
   % we'll expand this back to full space after filling it in
   dmq_map =  map_eps(ix_all(:,1), trip.Cannon.isArray);
   DMQ_eps = trip.eps(ix_all(:,1), trip.Cannon.isArray);
   SM1_eps = trip.eps(ix_all(:,2), trip.Cannon.isArray);
   SM2_eps = trip.eps(ix_all(:,3), trip.Cannon.isArray);

   DMQ_pvl = trip.pvl(ix_all(:,1), trip.Cannon.isArray);
   SM1_pvl = trip.pvl(ix_all(:,2), trip.Cannon.isArray);
   SM2_pvl = trip.pvl(ix_all(:,3), trip.Cannon.isArray);



   % append query interaction information: --- TAG ENTIRE ROWS ---------------------------------
   Qneg = Csv2Cell('DMQ_epsilon_neg.txt');
   neg_rows = ismember(trip.Cannon.Orf(ix_all(:,1)), Qneg);
   dmq_map(neg_rows,:) = cellfun(@appdQn, dmq_map(neg_rows,:), 'UniformOutput', false);

   Qpos = Csv2Cell('DMQ_epsilon_pos.txt');
   pos_rows = ismember(trip.Cannon.Orf(ix_all(:,1)), Qpos);
   dmq_map(pos_rows,:) = cellfun(@appdQp, dmq_map(pos_rows,:), 'UniformOutput', false);

   Qnan = Csv2Cell('DMQ_epsilon_nan.txt');
   nan_rows = ismember(trip.Cannon.Orf(ix_all(:,1)), Qnan);
   dmq_map(nan_rows,:) = cellfun(@appdQq, dmq_map(nan_rows,:), 'UniformOutput', false);


   
   % append array interaction information -------------------------------------
   DMQ_interactions = ((DMQ_eps <= -thresh) | (DMQ_eps >= thresh)) & (DMQ_pvl < 0.05);
   SMQ_interactions = ((SM1_eps <= -thresh | SM1_eps >= thresh) & (SM1_pvl < 0.05)) | ...
                      ((SM2_eps <= -thresh | SM2_eps >= thresh) & (SM2_pvl < 0.05));

   % additional information for printing
   SMQ_pos = (SM1_eps >= thresh & SM1_pvl < 0.05) | ...
             (SM2_eps >= thresh & SM2_pvl < 0.05);
   SMQ_neg = (SM1_eps <= -thresh & SM1_pvl < 0.05) | ...
             (SM2_eps <= -thresh & SM2_pvl < 0.05);

   % novel means, no overlap with significant digen
   novel_int = DMQ_interactions & ~SMQ_interactions;
   dmq_map(novel_int) = cellfun(@appdNOV, dmq_map(novel_int), 'UniformOutput', false);

   % All of these consider array interactions only
   quant_int = DMQ_interactions & SMQ_interactions; 
   quant_int_mix = DMQ_interactions & SMQ_neg & SMQ_pos; 
   quant_int_pos = DMQ_interactions & SMQ_pos; 
   quant_int_pos = quant_int_pos & ~quant_int_mix;
   quant_int_neg = DMQ_interactions & SMQ_neg; 
   quant_int_neg = quant_int_neg & ~quant_int_mix;


   % annotate all modified by array type
   dmq_map(quant_int_pos) = cellfun(@appdAp, dmq_map(quant_int_pos), 'UniformOutput', false);
   dmq_map(quant_int_neg) = cellfun(@appdAn, dmq_map(quant_int_neg), 'UniformOutput', false);
   dmq_map(quant_int_mix) = cellfun(@appdApn, dmq_map(quant_int_mix), 'UniformOutput', false);


   % rename Q-Novel to just Q-
   % dmq_map(strcmp('Q-Novel', dmq_map)) = {'Q-'};
   
   % -------------------------------------------------------------------

   % now count them:
   tri_types = unique(dmq_map(:));
   for i=1:length(tri_types);
      fprintf('%s\t%d\n', tri_types{i}, sum(sum(strcmp(tri_types{i}, dmq_map))));
   end

   % transcribe set info back to eps space
   map_eps(ix_all(:,1),trip.Cannon.isArray) = dmq_map;

end


function str = prepD(str)
   str = ['D' str];
end
function str = prepT(str)
   str = ['T' str];
end

function str = appdNOV(str)
   str = [str 'Novel'];
end

function str = appdQn(str)
   str = [str 'Q-'];
end
function str = appdQp(str)
   str = [str 'Q+'];
end
function str = appdQq(str)
   str = [str 'Q?'];
end

function str = appdAp(str)
   str = [str 'A+'];
end
function str = appdAn(str)
   str = [str 'A-'];
end
function str = appdApn(str)
   str = [str 'A+-'];
end



