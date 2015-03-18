function [pos_handle, neg_handle] = compare_sga_structs(sga, labels, standard, std_label)
%function [pos_handle, neg_handle] = compare_sga_structs(sga, labels, standard, std_label)

    PLOTS = true;
    SINGLE_PLOT = true;

    if(PLOTS)
        if(SINGLE_PLOT)
            neg_handle = subplot(2,2,3);
            pos_handle = subplot(2,2,4);
        else
            neg_handle = figure();
            pos_handle = figure();
        end
    end
    CM = cbrewer('qual', 'Set1', length(sga));

    for i=1:length(sga)
        %remove nans (not found in std)
        % mask out any queries or arrays not in the standard
        sga{i}.Cannon.isQuery(~ismember(StripOrfs(sga{i}.Cannon.Orf), standard.orfs)) = false;
        sga{i}.Cannon.isArray(~ismember(StripOrfs(sga{i}.Cannon.Orf), standard.orfs)) = false;

        % perform the analysis
        PR_neg = help_curve_the_pr(sga{i}, standard, true);
        PR_pos = help_curve_the_pr(sga{i}, standard, false);

        % draw the results
        if(PLOTS)
            % background, draw this for each SGA in case they are different.
            % background is the same for pos and negative.
            colors{i} = CM(i,:);
            backgrounds(i) = compute_background_rate(standard, sga{i});

            if(SINGLE_PLOT)
                subplot(neg_handle);
            else
                figure(neg_handle);
            end
            semilogx(PR_neg(:,3), PR_neg(:,2), 'Color', colors{i}, 'LineWidth', 2);
            xx = get(gca, 'XLim');
            set(gca, 'XLim', [min(xx(2)-1, 10), xx(2)]);
            hold on;

            if(SINGLE_PLOT)
                subplot(pos_handle);
            else
                figure(pos_handle);
            end

            semilogx(PR_pos(:,3), PR_pos(:,2), 'Color', colors{i}, 'LineWidth', 2);
            xx = get(gca, 'XLim');
            set(gca, 'XLim', [min(xx(2)-1, 10), xx(2)]);
            hold on;
        end
    end

    if(PLOTS)
        if(SINGLE_PLOT)
            subplot(neg_handle);
        else
            figure(neg_handle);
        end
        legend(labels, 'Interpreter', 'None');
        % Plot the backgrounds afterward so they don't show up in the legend
        for i=1:length(sga)
            xx = get(gca, 'XLim');
            plot(xx, [backgrounds(i) backgrounds(i)], ...
                'LineWidth', 2, 'LineStyle', '--', 'Color', colors{i});
        end
        xlabel('Recall');
        ylabel('Precision');
        title(['Negatives - ' std_label], 'Interpreter', 'None');
        hold off

        if(SINGLE_PLOT)
            subplot(pos_handle);
        else
            figure(pos_handle);
        end
        legend(labels, 'Interpreter', 'None');
        for i=1:length(sga)
            xx = get(gca, 'XLim');
            plot(xx, [backgrounds(i) backgrounds(i)], ...
                'LineWidth', 2, 'LineStyle', '--', 'Color', colors{i});
        end
        xlabel('Recall');
        ylabel('Precision');
        title(['Positives - ' std_label], 'Interpreter', 'None');
        hold off
    end
end

function [PR] = help_curve_the_pr(sga, standard, negativesTF)

   Qname = StripOrfs(sga.Cannon.Orf(sga.Cannon.isQuery));
   Aname = StripOrfs(sga.Cannon.Orf(sga.Cannon.isArray));
   EPS = sga.eps(sga.Cannon.isQuery, sga.Cannon.isArray);
   PVL = sga.pvl(sga.Cannon.isQuery, sga.Cannon.isArray);



   if(negativesTF)
      [r, c] = find(EPS <= -0.08 & PVL < 0.05);
   else
      [r, c] = find(EPS >=  0.08 & PVL < 0.05);
   end
   v = EPS(sub2ind(size(EPS), r, c));
   
   [tmp, ix] = sort(abs(v), 'descend');
   r = r(ix);
   c = c(ix);
   v = v(ix);

   [PR, ix] = help_eval_ppi_pairs(standard, [Qname(r) Aname(c)]);
   PR = [v(ix) PR];
end
function [PR, ix] = help_eval_ppi_pairs(standard, pair_list)

    name_map = Hash(java.util.HashMap(), standard.orfs);
   
    rows = zeros(size(pair_list, 1),1);
    cols = zeros(size(pair_list, 1),1);
    % none of these should return empty...
    for i=1:size(pair_list, 1); 
        rows(i) = name_map.get(pair_list{i,1});
        cols(i) = name_map.get(pair_list{i,2});
    end 
    linear_ind = sub2ind(size(standard.matrix), rows, cols);
    linear_vals = standard.matrix(linear_ind);
    if(islogical(standard.matrix)) % logical std
        ix = 1:length(linear_vals); %use them all, pass up the stack
        FP = cumsum(~linear_vals);
        TP = cumsum(linear_vals);
    else % +1 -1 std
        ix = find(linear_vals ~= 0); 
        TP = cumsum(linear_vals(ix) == 1); 
        FP = cumsum(linear_vals(ix) ==-1);
    end
    PR = [TP ./ (TP + FP), TP];
end
function [back] = compute_background_rate(standard, sga)
    % computes the rate of TP / (TP+FP) for SCREENED pairs
    % Assumes all queries are screened by all arrays
    % Also assumes standard is symmetric

    Q = ismember(standard.orfs, StripOrfs(sga.Cannon.Orf(sga.Cannon.isQuery)));
    A = ismember(standard.orfs, StripOrfs(sga.Cannon.Orf(sga.Cannon.isArray)));
    TP = sum(sum(standard.matrix(Q,A)== 1));
    if islogical(standard.matrix)
        TN = sum(sum(~standard.matrix(Q,A)));
    else
        TN = sum(sum(standard.matrix(Q,A)==-1));
    end
    back = TP / (TP+TN);

    % ? what p
end
