function [sgadata] = batch_correction_wrapper(sgadata, plate_id_map, lfid)
%% Batch correction


   % Get array plate means (use all screens, including WT screens)
   all_arrplates = unique(sgadata.arrayplateids);
   arrplate_ind = zeros(max(all_arrplates),1);
   arrplate_ind(all_arrplates) = 1:length(all_arrplates);  

   width = 48;
   height = 32;

   array_means = zeros(length(all_arrplates),width*height)+NaN;

   log_printf(lfid, ['Getting arrayplate means...\n|' blanks(50) '|\n|']);
   for i = 1:length(all_arrplates)
       
       currplates = unique(sgadata.plateids(sgadata.arrayplateids == all_arrplates(i)));
       t = [];
       
       for j = 1:length(currplates)
           
           ind = plate_id_map{currplates(j)};
           iii = sub2ind([32,48], sgadata.rows(ind), sgadata.cols(ind));
           
           d = zeros(32,48);
           d(:,[1 2 47 48]) = NaN;
           d([1 2 31 32],:) = NaN;
           d(iii) = sgadata.filt_colsize(ind);
           
           t = [t,d(:)];
           
       end
       
       array_means(i,:) = nanmean(t,2);
       
       % Print progress
       print_progress(lfid, length(all_arrplates),i);
       
   end
   log_printf(lfid, '|\n');

       

   % Do batch normalization using LDA
   save_mats=struct;

   log_printf(lfid, ['Preparing for batch normalization...\n|' blanks(50) '|\n|']);
   for i=1:length(all_arrplates)
       
       curr_plates = unique(sgadata.plateids(sgadata.arrayplateids == all_arrplates(i)));
       
       curr_mat = zeros(length(curr_plates),width*height)+NaN;
       curr_ind_mat = zeros(length(curr_plates),width*height)+NaN;
       curr_batch = zeros(length(curr_plates),1);
       
       for j = 1:length(curr_plates) 
           
           ind = plate_id_map{curr_plates(j)};
           iii = sub2ind([32,48], sgadata.rows(ind), sgadata.cols(ind));
           
           d = zeros(32,48);
           d(:,[1 2 47 48]) = NaN;
           d([1 2 31 32],:) = NaN;
           d_map = zeros(32,48) + 1; % default reference to the 1st index
           
           d_map(iii) = ind;
           d(iii) = sgadata.filt_colsize(ind);
           
           curr_mat(j,:)=d(:);
           curr_ind_mat(j,:)=d_map(:);
           
           curr_batch(j) = unique(sgadata.batch(ind));
           
       end
       
       t = curr_mat - repmat(array_means(i,:),size(curr_mat,1),1);
              
       save_mats(i).mat = t;
       save_mats(i).mat_ind = curr_ind_mat;
       save_mats(i).batch = curr_batch;
       
       % Print progress
       print_progress(lfid, length(all_arrplates),i);
       
   end
   log_printf(lfid, '|\n');
       

   % Normalize out batch effect. Method: LDA (supervised) 
   sgadata.batchnorm_colsize = sgadata.filt_colsize;
   log_printf(lfid, ['Batch normalization...\n|', blanks(50), '|\n|']);
   sv_thresh = 0.1;

   for i = 1:length(all_arrplates)
       
       t = save_mats(i).mat;
       t(isnan(t)) = 0;
       % batches_by_plate (formerly curr_batch) is a list of 
       % batch labels for all unique plateids in this array position
       batches_by_plate = save_mats(i).batch; 
     
       % Check if some of the batches are too small -- make a larger orphan batch
       % start assembling new batches by combining small batches until they reach certain size
       % TOO SMALL < 3
       % BIG ENOUGH >= 8 (as defined by median size in FG30)

       unique_batches_this_plate = unique(batches_by_plate);
       batch_count = histc(batches_by_plate, unique_batches_this_plate);

       merge_with_batch = find(batch_count < 3, 1, 'first'); % a pointer
       for j=1:length(unique_batches_this_plate)
           if(batch_count(j) < 3 && merge_with_batch ~= j) % don't merge batches with themselves 
               % merge this batch
               batches_by_plate(batches_by_plate == unique_batches_this_plate(j)) = ...
                                 unique_batches_this_plate(merge_with_batch);
               % update our counts and move our merge pointer if this orphan batch is big enough
               batch_count(merge_with_batch) = batch_count(merge_with_batch) + batch_count(j);
               batch_count(j) = NaN;
               if(batch_count(merge_with_batch) >=8 ) % move the pointer
                   merge_with_batch = merge_with_batch + find(batch_count(merge_with_batch+1:end) < 3, 1, 'first');
               end
           end
       end 
     
       batch_effect = multi_class_lda(t,batches_by_plate,sv_thresh);

       sgadata.batchnorm_colsize(save_mats(i).mat_ind(:)) = ...
          sgadata.batchnorm_colsize(save_mats(i).mat_ind(:)) - batch_effect(:);
       
       % Print progress
       print_progress(lfid, length(all_arrplates),i);
       
   end
   log_printf(lfid, '|\n');

end
