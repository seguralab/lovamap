% Extract information about 1D-ridge pathways that start at the center of
% the domain and radiate toward the edge
 
% peak_inds:          indices of peak points
% pks_graph_lengths:  list of rows indicating [peak 1, peak 2, ridge length]
%                     of peak graph
% edge_ind:           indices of edge of void space
% voxels:             list of [x,y,z] coordinates of all gridpoints in domain
% domain:             [xmin, xmax, ymin, ymax, zmin, zmax]
 
function [center_peak, path_nodes, path_length, path_r1Ds, path_tortuosity, path_neck_data, path_door_data, path_r1D_doors] = ...
               pathsCenter(peak_inds, pks_graph_length, ridges1D, ...
               edge_ind, voxels, domain)
 
    % Create a graph object for peaks and 1D-ridges
    pathway_graph = graph(pks_graph_length(:, 1), pks_graph_length(:, 2), ...
                          pks_graph_length(:, 3), 'omitselfloops');
    all_path_pks = unique([pks_graph_length(:, 1); pks_graph_length(:, 2)]);
    
    % locate true edge peaks
    edge_pks_touching = fastIntersect(edge_ind, peak_inds, 'elements');
    % convert to peak index
    edge_pks_ind = arrayfun(@(x) binarySearch(peak_inds, edge_pks_touching(x)), ...
                                 1 : numel(edge_pks_touching));
    % locate peaks connected to edge 1D-ridges
    [edge_pks_r1D, ind] = sort(ridges1D.rp_keyind(ridges1D.edge_bool, 1));
    if ~isempty(edge_pks_r1D)
        while edge_pks_r1D(1) == -1
            % remove cases of edge ridges that aren't flanked by any peaks
            edge_pks_r1D(1) = [];
            ind(1) = [];
        end
    end
    r1D_key = (1 : ridges1D.num)';
    r1D_key = r1D_key(ridges1D.edge_bool);
    r1D_key = r1D_key(ind);
    edge_pks = unique([edge_pks_r1D(:); edge_pks_ind(:)]);
    
    % find center of domain
    x_center = ((domain(2) - domain(1)) / 2) + domain(1);
    y_center = ((domain(4) - domain(3)) / 2) + domain(3);
    z_center = ((domain(6) - domain(5)) / 2) + domain(5);
    
    % narrow down peaks list to a cube surrounding center point
    cube_size = ((domain(2) - domain(1)) / 5) / 2; % 1/5 domain size
    x_range = [x_center - cube_size, x_center + cube_size];
    y_range = [y_center - cube_size, y_center + cube_size];
    z_range = [z_center - cube_size, z_center + cube_size];
    peaks_coord    = voxels(peak_inds, :); % already sorted
    peaks_cube_x   = peaks_coord(:, 1) >= x_range(1) & peaks_coord(:, 1) <= x_range(2);
    peaks_cube_y   = peaks_coord(:, 2) >= y_range(1) & peaks_coord(:, 2) <= y_range(2);
    peaks_cube_z   = peaks_coord(:, 3) >= z_range(1) & peaks_coord(:, 3) <= z_range(2);
    peaks_cube_ind = find(peaks_cube_x & peaks_cube_y & peaks_cube_z);
    % determine peak that is closest to center point of domain
    if ~isempty(peaks_cube_ind)
        [~, min_val] = min(arrayfun(@(x) norm(peaks_coord(peaks_cube_ind(x), :) - ...
                           [x_center, y_center, z_center]), 1 : length(peaks_cube_ind)));
        center_peak  = peaks_cube_ind(min_val);
    else
        [~, min_val] = min(arrayfun(@(x) norm(peaks_coord(x, :) - ...
                           [x_center, y_center, z_center]), 1 : length(peaks_coord)));
        center_peak  = min_val;
    end
    
    % Create key for pathway_graph edges and ridges1D
    [~, key] = sortrows(ridges1D.rp_keyind); % sort to match Edge table in pathway_graph
    r1D_pks2_sort = false(ridges1D.num, 1);
    r1D_pks2_sort(ridges1D.pks2) = true;
    r1D_pks2_sort = r1D_pks2_sort(key);
    key = key(r1D_pks2_sort);
    
    path_nodes  = cell(length(edge_pks) * 3, 1); % remove extra cells after
    path_length = zeros(length(edge_pks) * 3, 1);
    path_r1Ds  = cell(length(edge_pks) * 3, 1);
    path_tortuosity_lin = zeros(length(edge_pks) * 3, 1);
    path_tortuosity_vol = zeros(length(edge_pks) * 3, 1);
    path_tortuosity = struct;
    path_trkr = 1;
    for i = 1 : length(edge_pks)
        if binarySearch(all_path_pks, edge_pks(i)) > 0
            [a_node, b_length, c_edge] = shortestpath(pathway_graph, center_peak, edge_pks(i));
            if ~isempty(a_node)
                % Convert graph edge to ridge1D numbers
                c_edge = key(c_edge);
                % Add ridges1D that extend from edge_pk to edge of void space
                these_r1D = fastIntersect(edge_pks_r1D, edge_pks(i), 'indices');
                if ~isempty(these_r1D)
                    these_r1D = r1D_key(these_r1D);
                    % each additional ridge1D that extends to edge of void is a new path
                    for j = these_r1D(:)'
                        path_nodes{path_trkr}  = a_node;
                        path_length(path_trkr) = b_length;
                        path_r1Ds{path_trkr}  = c_edge;
                        % add edge path to form unique paths
                        path_r1Ds{path_trkr} = sort([path_r1Ds{path_trkr}(:); j]);
                        % add length
                        path_length(path_trkr) = path_length(path_trkr) + ridges1D.lengths(j);
                        % locate a ridge point at the edge of void space (end of path)
                        r1D_edgeind = fastIntersect(ridges1D.indices{j}, edge_ind, 'elements');
                        if ~isempty(r1D_edgeind)
                            % choose point that intersects with edge of void space
                            this_r1Dind = r1D_edgeind(1);
                        else
                            % if no indices intersect with edge of void space, scan for
                            % point that is physically farthest from edge_pk
                            phys_dist = arrayfun(@(x) norm(voxels(x, :) - voxels(peak_inds(edge_pks(i)), :)), ...
                                                      ridges1D.indices{j});
                            [~, max_ind] = max(phys_dist);
                            this_r1Dind = ridges1D.indices{j}(max_ind);
                        end
                        % Compute linear distance from edge to center for tortuosity
                        % measurement
                        lin_length = norm(voxels(this_r1Dind, :) - voxels(peak_inds(center_peak), :));
                        path_tortuosity_lin(path_trkr) = path_length(path_trkr) / lin_length;
                        % Compute volume of convex hull for another tortuosity measurement
                        comb_inds = [];
                        for k = 1 : length(path_r1Ds{path_trkr})
                            comb_inds = [comb_inds; ridges1D.indices{path_r1Ds{path_trkr}(k)}]; 
                        end
                        [~, conv_vol] = convhull(voxels(comb_inds, :));
                        path_tortuosity_vol(path_trkr) = conv_vol;
                        % update path number counter
                        path_trkr = path_trkr + 1;
                    end
                else
                    path_nodes{path_trkr}  = a_node;
                    path_length(path_trkr) = b_length;
                    path_r1Ds{path_trkr}  = c_edge;
                    % Compute linear distance from edge to center for tortuosity
                    % measurement
                    lin_length = norm(voxels(peak_inds(edge_pks(i)), :) - voxels(peak_inds(center_peak), :));
                    path_tortuosity_lin(path_trkr) = path_length(path_trkr) / lin_length;
                    % Compute volume of convex hull for another tortuosity measurement
                    comb_inds = [];
                    for j = 1 : length(path_r1Ds{path_trkr})
                        comb_inds = [comb_inds; ridges1D.indices{path_r1Ds{path_trkr}(j)}]; 
                    end
                    [~, conv_vol] = convhull(voxels(comb_inds, :));
                    path_tortuosity_vol(path_trkr) = conv_vol;
                    % update path number counter
                    path_trkr = path_trkr + 1;
                end
            end
        end
    end
    
    % Remove empty entries
    path_nodes = path_nodes(~cellfun('isempty', path_nodes));
    path_length = path_length(path_length ~= 0);
    path_r1Ds = path_r1Ds(~cellfun('isempty', path_r1Ds));
    path_tortuosity_lin = path_tortuosity_lin(path_tortuosity_lin ~= 0);
    path_tortuosity_vol = path_tortuosity_vol(path_tortuosity_vol ~= 0);

    % For each path, store the bottleneck widths (diameter) along 1D-ridges
    path_neck_data     = cell(length(edge_pks), 1);
    path_door_data     = cell(length(edge_pks), 1);
    path_r1D_doors = cell(length(edge_pks), 1);
    for i = 1 : numel(path_r1Ds)
        necks = zeros(numel(path_r1Ds{i}), 1);
        % make the following code more concise

        door_bool = ~fastIntersect(path_r1Ds{i}, ridges1D.connected, 'bool vec'); % finding 'doors' (i.e., the min on 1D ridges) that only exist between pores (not within pores)
        path_r1D_doors{i} = path_r1Ds{i}(door_bool);
        doorss = zeros(numel(path_r1D_doors{i}), 1);
        for j = 1 : numel(path_r1Ds{i})
            necks(j) = ridges1D.doors{path_r1Ds{i}(j)}.radius * 2;
        end
        for k = 1 : numel(path_r1D_doors{i})
            doorss(k) = ridges1D.doors{path_r1D_doors{i}(k)}.radius * 2;
        end
        path_neck_data{i} = necks;
        path_door_data{i} = doorss;
    end

    % save
    path_tortuosity.linear = path_tortuosity_lin;
    path_tortuosity.volume = path_tortuosity_vol;
    
    
%     % Starting at center peak, traverse peak graph outward
%     % store unique pathways that start at center point
%     pathways = mat2cell(peak_list_ind * ones(length(peak_inds), 1), ones(length(peak_inds), 1)); % over-allot and remove later
%     % find 1st layer of peaks connected to center peak
%     pathways_init = [find(adjacency(peak_list_ind, :)), find(adjacency(:, peak_list_ind))'];
%     % for 2nd layer of peaks, check which connecting peaks move outward
%     % Initialize eigenvectors
%     e_mat = [[1, 0, 0];
%              [-1, 0, 0];
%              [0, 1, 0];
%              [0, -1, 0];
%              [0, 0, 1];
%              [0, 0, -1]];
%     % each row in pathways_q is the queue to a different pathway
%     pathways_q = pathways_init(:);
%     which_pathway = 1;
%     for i = 1 : length(pathways_init)
%         pathways{which_pathway} = [pathways{which_pathway}, pathways_init(i)];
%         keep_going = 1;
%         while pathways_q(i) ~= 0
%             check_pks = [find(adjacency(pathways_init(i), :)), find(adjacency(:, pathways_init(i)))'];
%             if ~isempty(check_pks)
%                 for j = check_pks(:)'
%                     [~, closest_wall] = max(arrayfun(@(x) dot(voxels(peak_inds(j), :), e_mat(x, :)), ...
%                                             1 : 6));
%                     % dot product to closest e_vector of closest surface
%                     direc_vec = dot(voxels(peak_inds(j), :), e_mat(closest_wall, :));
%                     if direc_vec > 0
%                         % peak is not heading back toward center
%                         % add peak to pathways list
%                         pathways{which_pathway} = [pathways{which_pathway}, j];
%                         pathways_q(i) = j;
%                     end
%                 end
%             else
%                 pathways_q(i) = 0;
%             end
%              
% %             switch closest_wall
% %                 case 1
% %                     dist = domain(2) - voxels(peak_inds(j), 1);
% %                 case 2
% %                     dist = voxels(peak_inds(j), 1) - domain(1);
% %                 case 3
% %                     dist = domain(4) - voxels(peak_inds(j), 2);
% %                 case 4
% %                     dist = voxels(peak_inds(j), 2) - domain(3);
% %                 case 5
% %                     dist = domain(6) - voxels(peak_inds(j), 3);
% %                 case 6
% %                     dist = voxels(peak_inds(j), 3) - domain(5);
% %                 otherwise
% %                     error('Are you for sure in 3-D?')
% %             end
% 
%         end
%     end
%     
%     
%     path_cntr  = 1;
%     while ~isempty(pathways_q)
%         mrkr_pk = pathways_q(end);
%         pathways{path_cntr} = [pathways{path_cntr}, mrkr_pk];
%         check_pks = [find(adjacency(mrkr_pk, :)), find(adjacency(:, mrkr_pk))];
%         
%         
%     end
    
end
