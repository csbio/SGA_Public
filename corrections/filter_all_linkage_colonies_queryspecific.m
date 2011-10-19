%%
% FILTER_ALL_LINKAGE_COLONIES_QUERYSPECIFIC - generates the list of colonies affected by linkage
%
% Inputs:
%   sgadata - structure containing all the colony size data
%   linkagefile - the name of the file containing the coordinates of the linkage windows
%
% Outputs:
%   all_linkage_cols - indices of the colonies affected by linkage
%
% Authors: Chad Myers (cmyers@cs.umn.edu), 
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca),
%          Benjamin VanderSluis (bvander@cs.umn.edu)
%
% Last revision: 2011-09-19
%
%%

function all_linkage_cols = filter_all_linkage_colonies_queryspecific(sgadata, linkagefile)

    % Print the name and path of this script
    p = mfilename('fullpath');
    fprintf('\nLinkage filter script:\n\t%s\n\n',p);
    
    % Load chromosomal coordinates

    %fprintf('Using standard chromosome coordinates\n');
    %chromdata = importdata('chrom_coordinates.txt');

    fprintf('Using TSA chromosome coordinates\n');
    chromdata = importdata('chrom_coordinates_wTS_110907.txt'); %ek 110907

    orfs_coords = chromdata.textdata;
    chrom_coords = [chromdata.data(:,1),chromdata.data(:,2),chromdata.data(:,3)];

    all_linkage_cols = []; 
    linkage_dist = 200e3;   % default linkage distance
    
    % Load query-specific linkage
    % Linkage file will now accept NANs
    [lnkg.textdata, lnkg.data(:,1), lnkg.data(:,2)] = textread(linkagefile,'%s %d %d');
    lnkg.data = sort(lnkg.data, 2, 'ascend'); % some regions are defined high to low, reverse them
    linkage_reg = zeros(length(sgadata.orfnames),4) - 1; %BJV 4 cols
    [t,ind1,ind2] = intersect(sgadata.orfnames, lnkg.textdata);
    linkage_reg(ind1,[1,2]) = lnkg.data(ind2,:);

    % If no single query-specific linkage info is available, assume 200kb
    all_querys = unique(sgadata.querys); % BJV this is just  1:max(sgadata.querys) ie 1:1713
    ind = find(linkage_reg(all_querys,1) < 0);
    fprintf('No pre-defined linkage info for these queries (200 kb have to be assumed):\n\tDouble Queries will be parsed further\n');
    fprintf('%s\n', sgadata.orfnames{all_querys(ind)});
    
    fprintf('\nNo coordinate or linkage info for these queries:\n');
    for i = 1:length(ind)
        orf_string = sgadata.orfnames{all_querys(ind(i))};
        split_genes_char = findstr(orf_string, '+'); % may be double
        assert(length(split_genes_char) <= 1); % assumption: [0,1] +'s
        if(isempty(split_genes_char)) % single gene
            orf1 = orf_string;
            orf2 = '';

            split_ann_char = findstr(orf1, '_');
				if(~isempty(split_ann_char))
					orf1 = orf1(1:split_ann_char-1); % remove annotation
				end
            
            c1 = strmatch(orf1, orfs_coords, 'exact');
            l1 = strmatch(orf1, lnkg.textdata, 'exact');
            c2 = []; 
            l2 = []; 

        else % maybe a double gene
            orf1 = orf_string(1:split_genes_char-1);
            orf2 = orf_string(split_genes_char+1:end);

            split_ann_char = findstr(orf1, '_');
				if(~isempty(split_ann_char))
					orf1 = orf1(1:split_ann_char-1); % remove annotation
				end

            split_ann_char = findstr(orf2, '_');
				if(~isempty(split_ann_char))
					orf2 = orf2(1:split_ann_char-1); % remove annotation
				end

            c1 = strmatch(orf1, orfs_coords, 'exact');
            l1 = strmatch(orf1, lnkg.textdata, 'exact');
            c2 = strmatch(orf2, orfs_coords, 'exact');
            l2 = strmatch(orf2, lnkg.textdata, 'exact');
        end
       
        % Handle the first "ORF"
        if(~isempty(l1)) % single ORF 1 found in file
            linkage_reg(all_querys(ind(i)),[1,2]) = lnkg.data(l1,:);
        elseif(~isempty(c1)) % insert default for this ORF
            linkage_reg(all_querys(ind(i)),[1,2]) = [max(chrom_coords(c1,2)-linkage_dist,1),chrom_coords(c1,3)+linkage_dist];
        else
            fprintf('E1 no linkage or coord info for %s in %s\n', orf1, orf_string); % first str empty
        end

        % Handle the second "ORF"
        if(~isempty(l2)) % single ORF 2 found in file
            linkage_reg(all_querys(ind(i)),[3,4]) = lnkg.data(l2,:);
        elseif(~isempty(c2)) % insert default for this ORF
            linkage_reg(all_querys(ind(i)),[3,4]) = [max(chrom_coords(c2,2)-linkage_dist,1),chrom_coords(c2,3)+linkage_dist];
        else
				if(~isempty(orf2)) % Don't complain about no info for an empty orf (see line 53: orf2 = '')
            	fprintf('E2 no linkage or coord info for %s in %s\n', orf2, orf_string);
				end
        end

    end

    % URA3 linkage (WT-specific)
    t = get_linked_ORFs('YEL021W', orfs_coords, linkage_dist, chrom_coords);
    [int,ia,ib] = intersect(t, sgadata.orfnames);
    ura3_linkage_cols = find(ismember(sgadata.arrays,ib) & ismember(sgadata.querys,strmatch('undefined',sgadata.orfnames)));
    all_linkage_cols = [all_linkage_cols; ura3_linkage_cols];
    
    % HO linkage (Triple WT-specific) for now we'll pull this from ALL WT queries
    % This should be handled from code above
    %{
    t = get_linked_ORFs('YDL227C', orfs_coords, linkage_dist, chrom_coords);
    [int,ia,ib] = intersect(t, sgadata.orfnames);
    ho_linkage_cols = find(ismember(sgadata.arrays,ib) & ismember(sgadata.querys,strmatch('undefined+YDL227C',sgadata.orfnames)));
    all_linkage_cols = [all_linkage_cols; ho_linkage_cols];
    %}

    % LYP1 linkage
    t = get_linked_ORFs('YNL268W', orfs_coords, linkage_dist, chrom_coords);
    [int,ia,ib] = intersect(t, sgadata.orfnames);
    lyp1_linkage_cols = find(ismember(sgadata.arrays,ib));
    all_linkage_cols = [all_linkage_cols; lyp1_linkage_cols];
    
    % CAN1 linkage
    t = get_linked_ORFs('YEL063C', orfs_coords, linkage_dist, chrom_coords);
    [int,ia,ib] = intersect(t, sgadata.orfnames);
    can1_linkage_cols = find(ismember(sgadata.arrays,ib));
    all_linkage_cols = [all_linkage_cols; can1_linkage_cols];

    exp_linkage_cols = [];
    
    fprintf(['Mapping query-specific linkage...\n|', blanks(50), '|\n|']);
    y = 0;
    
    for i = 1:length(all_querys)
        
        %% First Linkage
        lnkg_reg = linkage_reg(all_querys(i),[1,2]); % first linkage; all_querys(i) == i BJV
        indtt = [];
        t = {};
        if(lnkg_reg(1) ~= -1)
            curr_chrom = findstr(sgadata.orfnames{all_querys(i)}(2), 'ABCDEFGHIJKLMNOP');
            if isempty(curr_chrom)
                if ~strcmp('undefined',sgadata.orfnames{all_querys(i)})
                    fprintf('E1 Could not get a chromosome # for this orf: %s\n', sgadata.orfnames{all_querys(i)});
                end
                continue;
            end
            
            % Get query-specific linked ORFs (= a set of ORFs that start or end within the linkage region)
            indt = find(chrom_coords(:,1) == curr_chrom);
    
            indtt = find((min(chrom_coords(indt,2:3),[],2) >= min(lnkg_reg) & min(chrom_coords(indt,2:3),[],2) < max(lnkg_reg)) | ...
            (max(chrom_coords(indt,2:3),[],2) >= min(lnkg_reg) & max(chrom_coords(indt,2:3),[],2) < max(lnkg_reg)));

            indtt = unique(indtt); 
            t = orfs_coords(indt(indtt));
        end
	
        %% Second Linkage
        lnkg_reg = linkage_reg(all_querys(i),[3,4]); % second linkage; all_querys(i) == i BJV
        indtt_double = [];
        t_d = {};
        if(lnkg_reg(1) ~= -1) % skip this step for single_deletion queries BJV
            % ASSUMPTION: if we get here then we've already successfully split this string
            whole_string = sgadata.orfnames{all_querys(i)};
            split_char = findstr('+', whole_string);
            curr_chrom = findstr(whole_string(split_char+2), 'ABCDEFGHIJKLMNOP');
            if isempty(curr_chrom)
                if ~strcmp('undefined',sgadata.orfnames{all_querys(i)})
                    fprintf('E2 Could not get a chromosome # for this orf: %s\n', sgadata.orfnames{all_querys(i)});
                end
                continue;
            end
    
            indt_double = find(chrom_coords(:,1) == curr_chrom);

            indtt_double = find((min(chrom_coords(indt_double,2:3),[],2) >= min(lnkg_reg) & min(chrom_coords(indt_double,2:3),[],2) < max(lnkg_reg)) | ...
            (max(chrom_coords(indt_double,2:3),[],2) >= min(lnkg_reg) & max(chrom_coords(indt_double,2:3),[],2) < max(lnkg_reg)));
    
            indtt_double = unique(indtt_double);
            t_d = orfs_coords(indt_double(indtt_double));
        end

        t = unique([t; t_d]); % merge the two orf lists
        
        [int,ia,ib] = intersect(t, sgadata.orfnames);
                
        curr_inds = find(ismember(sgadata.arrays,ib) & sgadata.querys == all_querys(i));
        exp_linkage_cols = [exp_linkage_cols; curr_inds];
        
        % Print progress
        x = fix(i * 50/length(all_querys));
        if x > y
            fprintf('*')
            y = x;
        end
    end
    fprintf('|\n');

    all_linkage_cols = [all_linkage_cols; exp_linkage_cols];
