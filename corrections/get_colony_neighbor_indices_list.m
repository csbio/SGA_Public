%%
% GET_COLONY_NEIGHBOR_INDICES_LIST - generates a list of indices to plate positions
%   corresponding to 3 colonies from the 3 neighboring strains
%
% Inputs:
%   none
%
% Outputs:
%   colony_neighbors_list - corrected colony sizes
%
% Authors: Chad Myers (cmyers@cs.umn.edu), Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
%
% Last revision: 2010-07-19
%
%%

function colony_neighbors_list = get_colony_neighbor_indices_list()

    % Print the name and path of this script
    p = mfilename('fullpath');
    fprintf('\n\tGet colony neighbor indices script:\n\t\t%s\n\n',p);

    tmp_inds = zeros(32,48);
    tmp_inds(:)=1:length(tmp_inds(:));

    colony_neighbors_list = zeros(1536,3)+1;


    for i = 3:30
        for j=3:46
            
            if mod(i,2) == 1 % odd rows
                if mod(j,2) == 1 % odd columns
                    colony_neighbors_list(tmp_inds(i,j),:) = [tmp_inds(i,j-1),tmp_inds(i-1,j),tmp_inds(i-1,j-1)];
                else % even columns
                    colony_neighbors_list(tmp_inds(i,j),:) = [tmp_inds(i,j+1),tmp_inds(i-1,j),tmp_inds(i-1,j+1)];
                end
            else % even rows
                if mod(j,2) == 1 % odd columns
                    colony_neighbors_list(tmp_inds(i,j),:) = [tmp_inds(i,j-1),tmp_inds(i+1,j),tmp_inds(i+1,j-1)];
                else % even columns
                    colony_neighbors_list(tmp_inds(i,j),:) = [tmp_inds(i,j+1),tmp_inds(i+1,j),tmp_inds(i+1,j+1)];
                end            
            end

        end
    end


