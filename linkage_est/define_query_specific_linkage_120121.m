% Input: 
%   dataset.queries     Cell array of query ORFs (with or without annot.)
%   dataset.arrays      Cell array of array ORFs (with or without annot.)
%   dataset.scores  Matrix containing the genetic interaction scores
%
% Output:
%   lnkg.orf            Cell array of query ORFs (as in input)
%   lnkg.coord_mean          Matrix (queries x 2) with the linkage boundaries

function [dataset_filt, dataset_corr, dataset_noncorr, lnkg] = define_query_specific_linkage_120121(dataset)


%% Load the data

[t,ix1] = sort_by_chr(dataset.queries);
[t,ix2] = sort_by_chr(dataset.arrays);
dataset.queries = dataset.queries(ix1);
dataset.arrays = dataset.arrays(ix2);
dataset.scores = dataset.scores(ix1,ix2);

%dataset.scores(isnan(dataset.scores)) = 0;

labels_row = dataset.queries;
labels_col = dataset.arrays;
data = dataset.scores;


chromosomes = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'};

screens = labels_row;
[labels_row, annotation] = strtok(labels_row,'_');

arrays = labels_col;
[labels_col, annotation] = strtok(labels_col,'_');

load chr_length_110207;

%load orf_coordinates_110207;
load orf_coordinates_150617;
coord = orf_coord; clear orf_coord;

%% Get the chromosomal coordinates of queries and arrays
inds1 = [1 : length(labels_row)]';
inds2 = multistrmatch(labels_row, coord.orf,1,1,1);
ii = find(inds2 == 0);
inds1(ii) = []; 
inds2(ii) = [];

labels_row = labels_row(inds1);
data = data(inds1,:);
coord_row = [coord.start(inds2) coord.end(inds2)];

screens = screens(inds1);

inds1 = [1 : length(labels_col)]';
inds2 = multistrmatch(labels_col, coord.orf,1,1,1);
ii = find(inds2 == 0);
inds1(ii) = [];
inds2(ii) = [];

labels_col = labels_col(inds1);
data = data(:,inds1);
coord_col = [coord.start(inds2) coord.end(inds2)];

arrays = arrays(inds1);

lnkg.orf = screens;
lnkg.coord_mean = zeros(length(screens),2)+NaN;
%lnkg.coord_median = zeros(length(labels_row),2)+NaN;

dataset_filt.queries = screens;
dataset_filt.arrays = arrays;
dataset_filt.scores = data;

dataset_corr.queries = screens;
dataset_corr.arrays = arrays;
dataset_corr.scores = zeros(length(screens), length(arrays))+NaN;

dataset_noncorr.queries = screens;
dataset_noncorr.arrays = arrays;
dataset_noncorr.scores = zeros(length(screens), length(arrays))+NaN;


%%

redfactor = 100;
for ic = 1:16
    
    iq = strmatch(['Y' chromosomes(ic)], labels_row);
    ia = strmatch(['Y' chromosomes(ic)], labels_col);
    
    all_data(ic).labels_row = labels_row(iq);
    all_data(ic).labels_col = labels_col(ia);
    all_data(ic).coord_row = round(coord_row(iq,:)/redfactor);
    all_data(ic).coord_col = round(coord_col(ia,:)/redfactor);
    all_data(ic).scores = data(iq,ia);
    all_data(ic).coord_col_length = round(coord_col(ia,2)/redfactor)-round(coord_col(ia,1)/redfactor)+1;
    all_data(ic).coord_col_mid = all_data(ic).coord_col(:,1) + round((all_data(ic).coord_col(:,2)-all_data(ic).coord_col(:,1))/2);
    all_data(ic).coord_row_mid = all_data(ic).coord_row(:,1) + round((all_data(ic).coord_row(:,2)-all_data(ic).coord_row(:,1))/2);
    
    linkage_profile = zeros(length(iq), round(chr_length.length(ic)/redfactor));
    linkage_profile_unscaled = zeros(length(iq), round(chr_length.length(ic)/redfactor));
    
    % The size of the window being averaged (in bp)
    smooth_window2 = round(60000/redfactor);
    
    % Quite complicated way to speed up the smooth process:
    % Create a matrix of arrays x chromosomal positions where each element
    % indicates how much of a particular array is covered by a window
    % centered on that chromosomal position
    
    coverage.bp = zeros(length(ia),round(chr_length.length(ic)/redfactor));
    
    
    L = repmat(all_data(ic).coord_col_length, 1, round(chr_length.length(ic)/redfactor));
    S = repmat(all_data(ic).coord_col(:,1), 1, round(chr_length.length(ic)/redfactor));
    E = repmat(all_data(ic).coord_col(:,2), 1, round(chr_length.length(ic)/redfactor));
    X = repmat(1:round(chr_length.length(ic)/redfactor), length(ia),1);
    M = repmat(all_data(ic).coord_col_mid, 1, round(chr_length.length(ic)/redfactor));
    
    W2 = smooth_window2/2;

    fprintf(['\nChr ' num2str(ic) ' - Function in progress...\n|' blanks(50) '|\n|']);
    
    tot_iq = length(iq);
    for q = 1 : tot_iq
        
        ix = find(all_data(ic).scores(q,:) > 40); all_data(ic).scores(q,ix) = 40;   % eliminate the most positive scores to reduce the noise
        
        D = all_data(ic).coord_row_mid(q) - X;        
                
        coverage.bp = L - max(0, E-(X+W2)) - max(0, (X-W2)-S);
        coverage.bp(coverage.bp<0) = 0;
        
         % Ignore the most positive score in a window
         scores = repmat(all_data(ic).scores(q,:)',1,round(chr_length.length(ic)/redfactor)); 
         scores(coverage.bp==0) = NaN;
         %m = repmat(max(scores,[],1),length(ia),1);
         %mn = min(scores,[],1);
         mx = repmat(max(scores,[],1),length(ia),1);
         
         scores(scores==mx) = NaN;
        
%         % Give more weight to the numbers in the middle of the window
        all_windows_mid = repmat(1:round(chr_length.length(ic)/redfactor),length(ia),1);
        all_arrays_mid = repmat(all_data(ic).coord_col_mid,1,round(chr_length.length(ic)/redfactor));
        t1 = 1-abs(all_windows_mid-all_arrays_mid)./W2;
        
        % Give more weight to the numbers closer to the query
%         all_arrays_mid = repmat(all_data(ic).coord_col_mid,1,round(chr_length.length(ic)/redfactor));
%         d = abs(all_arrays_mid - all_data(ic).coord_row_mid(q));
%         d1 = d;
%         d1(coverage.bp==0) = NaN;
%         mn = repmat(min(d1,[],1), length(ia),1);
%         mx = repmat(max(d1,[],1), length(ia),1);
%         
%         r = [1 0.1];
%         
%         t1 = r(2) + (r(1)-r(2)) .* (d1 - mx) ./ (mn - mx);
        
        
        % Limit the windows next to the query to the smallest possible
        % width
        inds = find(abs(D) < W2);
        coverage.bp(inds) = L(inds) - max(0, E(inds)-(X(inds)+abs(D(inds)))) - max(0, (X(inds)-abs(D(inds)))-S(inds));

        inds = find(D<0 & D>-W2 & M>X);
        coverage.bp(inds) = 0;
        
        inds = find(D>0 & D<W2 & M<X);
        coverage.bp(inds) = 0;
        
        
        coverage.bp(coverage.bp<0) = 0;
        coverage.bp(isnan(scores)) = 0;
       
        
        coverage.bp_binary = coverage.bp;
        coverage.bp_binary(coverage.bp_binary>0)=1;

        
        coverage.bp_weighted = coverage.bp .* t1;
        
              
        % Calculating the mean        
        
        linkage_profile(q,:) = nansum(coverage.bp_weighted .* scores,1)./nansum(coverage.bp_weighted,1);
        linkage_profile_unscaled(q,:) = nansum(coverage.bp .* scores,1)./nansum(coverage.bp,1);
               
        % Calculating the linkage window boundaries
        boundary_left2 = max(find(linkage_profile(q,1:all_data(ic).coord_row(q,1)) >= 0));

        if isempty(boundary_left2)
            boundary_left2 = 1;
        end
        
        boundary_right2 = all_data(ic).coord_row(q,2) + min(find(linkage_profile(q, all_data(ic).coord_row(q,2):end)>=0));

        if isempty(boundary_right2)
            boundary_right2 = round(chr_length.length(ic)/redfactor);
        end
        
    %% Figure plotting. If uncommented, will create a plot for each query. Good to try out the first time you run a linkage analysis on a new dataset.
%         figure(iq(q));
%          %subplot(2,1,1);
%         hold on;
%         for a = 1 : length(ia)
%             patch([all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,2) all_data(ic).coord_col(a,2)],[0 all_data(ic).scores(q,a) all_data(ic).scores(q,a) 0],'k');
%         end
%         plot(1:round(chr_length.length(ic)/redfactor), linkage_profile(q,:),'m-','LineWidth',3);
%        % plot(1:round(chr_length.length(ic)/redfactor), linkage_profile_unscaled(q,:),'y-','LineWidth',2);
%         %plot(1:chr_length.length(ic), profile(q).median,'c-','LineWidth',3);
%         
%         plot([boundary_left2 boundary_left2], [min(linkage_profile(q,:)) max(linkage_profile(q,:))],'g-','LineWidth',3); 
%         plot([boundary_right2 boundary_right2], [min(linkage_profile(q,:)) max(linkage_profile(q,:))],'g-','LineWidth',3);
%         
% 
%         plot([chr_length.length(ic)-smooth_window2 chr_length.length(ic)],[-1 -1], 'k-','LineWidth',2);
%         %plot([boundary_right-smooth_window2/2 boundary_right+smooth_window2/2],[-100 -100],'k-','LineWidth',2);
%         grid on;
%         set(gca,'XLim',[1 round(chr_length.length(ic)/redfactor)]);
%         
%         title(screens(iq(q)),'Interpreter','none');
%         linkage_profile(q, 1:boundary_left2-1) = 0;
%         linkage_profile(q, boundary_right2+1:end) = 0;
%         
%         m = all_data(ic).coord_row(q,1) + (all_data(ic).coord_row(q,2)-all_data(ic).coord_row(q,1))/2;
%         plot([m-smooth_window2/2 m+smooth_window2/2],[1 1],'k-');
%         plot(m,1,'ko','MarkerFaceColor','k');
      
                
    %% Define linkage boundaries
    
        bl2 = boundary_left2*redfactor;
        br2 = boundary_right2*redfactor;
        
        if bl2 == redfactor
            bl2 = 1;
        end
        
        if br2 > chr_length.length(ic)
            br2 = chr_length.length(ic);
        end
        
        lnkg.coord_mean(iq(q),:) = [bl2 br2];

    %% Identify the arrays that are in the linkage window
        
        if isempty(find(coverage.bp_binary > 0))
            i1 = 1;
            i2 = 1;
        else
        
            i1 = []; step = 0;
            while isempty(i1)
                i1 = min(find(coverage.bp_binary(:,boundary_left2+step) > 0));
                step = step + 100;
            end

            i2 = []; step = 0;
            while isempty(i2) & boundary_right2+step <= size(coverage.bp_binary,2)
                i2 = max(find(coverage.bp_binary(:,boundary_right2+step) > 0));
                step = step - 100;
            end
            
        end
        
        dataset_filt.scores(iq(q),ia(i1:i2)) = NaN;
        
    %% Correct within-linkage scores by the estimated profile
        
%         subplot(2,1,2);
%         hold on;
%         for a = 1 : i1-1
%             patch([all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,2) all_data(ic).coord_col(a,2)],[0 all_data(ic).scores(q,a) all_data(ic).scores(q,a) 0],'k');
%         end
        for a = i1 : i2
            
            y = all_data(ic).scores(q,a) - min(linkage_profile(q,all_data(ic).coord_col(a,1):all_data(ic).coord_col(a,2)));
            dataset_corr.scores(iq(q), ia(a)) = y;
            dataset_noncorr.scores(iq(q),ia(a)) = all_data(ic).scores(q,a);
            
%             patch([all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,2) all_data(ic).coord_col(a,2)],[0 y y 0],'k');

        end
%         for a = i2+1 : length(ia)
%             patch([all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,1) all_data(ic).coord_col(a,2) all_data(ic).coord_col(a,2)],[0 all_data(ic).scores(q,a) all_data(ic).scores(q,a) 0],'k');
%         end
%         
%         plot([boundary_left2 boundary_left2], [min(linkage_profile(q,:)) max(linkage_profile(q,:))],'r-','LineWidth',3); 
%         plot([boundary_right2 boundary_right2], [min(linkage_profile(q,:)) max(linkage_profile(q,:))],'r-','LineWidth',3);
%         
%         grid on;
%         set(gca,'XLim',[1 round(chr_length.length(ic)/redfactor)]);
        
%         m = all_data(ic).coord_row(q,1) + (all_data(ic).coord_row(q,2)-all_data(ic).coord_row(q,1))/2;
        %plot([m-smooth_window2/2 m+smooth_window2/2],[50 50],'k-');
        %plot(m,50,'ko','MarkerFaceColor','k');
        
%%
        print_progress(1, length(iq), q);
        
        
        
    end
    
    
%     query_linkage_profile.query = labels_row(iq);
%     query_linkage_profile.profile_sparse = sparse(linkage_profile);
%     clear linkage_profile;
% 
%     save(['Laboratory/SGA/Recombination_map/11-09-23/query_linkage_profile_' num2str(ic) '_sparse.mat'],'query_linkage_profile','-v7.3');
%     clear query_linkage_profile;
    
    fprintf('|\n');
end

    
