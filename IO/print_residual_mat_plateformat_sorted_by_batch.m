function print_residual_mat_plateformat_sorted_by_batch(sgadata,field,filename,query_map, plate_id_map)

avg_mats=struct;
avgt=[];
all_arrplates = unique(sgadata.arrayplateids);
width = 48;
height = 32;

for i=1:length(all_arrplates)
    
    curr_plates = unique(sgadata.plateids(sgadata.arrayplateids == all_arrplates(i)));
    
    curr_mat = zeros(length(curr_plates),width*height)+NaN;
    curr_ind_mat = zeros(length(curr_plates),width*height)+NaN;
    curr_batch = [];
    t=[];

    
    for j=1:length(curr_plates)        
        
        % Get a matrix of colony sizes (d) and colony indices (d_map) for the plate
        ind = plate_id_map{curr_plates(j)};
        d = zeros(32,48);
        d(:,[1,2,47,48]) = NaN;
        d([1,2,31,32],:) = NaN;
        
        d_map = zeros(32,48)+1; % Default reference to 1st index
        
        iii = sub2ind([32,48], sgadata.rows(ind), sgadata.cols(ind));
        d_map(iii) = ind;
        d(iii) = sgadata.(field)(ind);
        
        curr_mat(j,:)=d(:);
        curr_ind_mat(j,:)=d_map(:);
        
        t=d;
        %if(mod(j,10)==0)  fprintf('Finished %d of %d...\n',j,length(curr_plates)); end        
        %t=[t;zeros(2,size(t,2))];
        avg_mats(i,j).mat=t;
    end
    
    fprintf('Finished plate %d...\n', i);
    
end
    

for i=1:size(avg_mats,1)
    
   sum=zeros(size(avg_mats(1,1).mat));
   num=zeros(size(avg_mats(1,1).mat));
   
    for j=1:size(avg_mats,2),
        ind = find(~isnan(avg_mats(i,j).mat));
        sum(ind) = sum(ind) + avg_mats(i,j).mat(ind);
        num(ind) = num(ind)+1;
    end
   
    comb_avg(i).mat = sum./num;
end

for i=1:length(all_arrplates)
    
    curr_plates = unique(sgadata.plateids(sgadata.arrayplateids == all_arrplates(i)));
    lp = length(curr_plates);
    
    curr_mat = zeros(length(curr_plates),width*height)+NaN;
    curr_ind_mat = zeros(length(curr_plates),width*height)+NaN;
    curr_batch = [];
    t=zeros(34*lp, 48);
    save_mats(i).querys=zeros(1,34*lp)-1;
    save_mats(i).rows=zeros(1,34*lp)-1;
    save_mats(i).setids=zeros(1,34*lp)-1;
    
    c1 = 1;
    for j=1:lp          
        
        inds_p = find(sgadata.plateids == curr_plates(j));
        
        % Get a matrix of colony sizes (d) and colony indices (d_map) for the plate
        ind = plate_id_map{curr_plates(j)};
        d = zeros(32,48);
        d(:,[1,2,47,48]) = NaN;
        d([1,2,31,32],:) = NaN;
        
        d_map = zeros(32,48)+1; % Default reference to 1st index
        
        iii = sub2ind([32,48], sgadata.rows(ind), sgadata.cols(ind));
        d_map(iii) = ind;
        d(iii) = sgadata.(field)(ind);
        
        curr_mat(j,:)=d(:);
        curr_ind_mat(j,:)=d_map(:);

        d = d - comb_avg(i).mat;

        t(c1:c1+31,:)=d;
        

        if(mod(j,10)==0)  fprintf('Finished %d of %d...\n',j,lp); end
        save_mats(i).querys(c1:c1+31) = unique(sgadata.querys(inds_p));
        save_mats(i).rows(c1:c1+31) = [1:32];      
        save_mats(i).setids(c1:c1+31) = unique(sgadata.setids(inds_p));
        
        %t=[t;zeros(2,size(t,2))];
        
%         save_mats(i).querys(c1+33:c1+34)=-1;
%         save_mats(i).rows(c1+33:c1+34)=-1;
%         save_mats(i).setids(c1+33:c1+34)=-1;
        
        c1 = c1+34;
    end
    
    save_mats(i).mat =t;
    
    
end
    


xstrs = [];
ystrs=[];



curr_mat=[];
all_batches = unique(sgadata.batch);

count=1;
ystrs={};
for k = 1 : length(all_batches)
    inds = find(sgadata.batch == all_batches(k));
    
    all_querys = unique(sgadata.querys(inds));
    
    for i=1:length(all_querys)
        %check if there are multiple sets
        curr_sets = unique(sgadata.setids(intersect(query_map{all_querys(i)}, inds)));

        for k=1:length(curr_sets)
            for j=1:length(all_arrplates)
               if(~isempty(find(save_mats(j).querys == all_querys(i) & save_mats(j).setids == curr_sets(k))) & length(find(save_mats(j).querys == all_querys(i) & save_mats(j).setids==curr_sets(k))) <= 32)
                   curr_mat(((count-1)*32+1):((count-1)*32+32),((j-1)*48+1):((j-1)*48+48)) = save_mats(j).mat(save_mats(j).querys == all_querys(i) & save_mats(j).setids==curr_sets(k),:);
               else
                   curr_mat(((count-1)*32+1):((count-1)*32+32),((j-1)*48+1):((j-1)*48+48)) = zeros(height,width)+NaN;
               end
            end

            count = count+1;
            for j=1:32,
                ystrs = [ystrs; {['Query: ',sgadata.orfnames{all_querys(i)},' Set: ', num2str(curr_sets(k)), ' Row: ',num2str(j)]}];
            end

        end


    end
end

count=1;
xstrs = [];

for j=1:length(all_arrplates),
    for i=1:48,
        xstrs{count} = ['Array plate ',num2str(j),'; Column ',num2str(i)];
        count=count+1;
    end
end


print_pcl_file(curr_mat,ystrs,ystrs,xstrs,filename),
