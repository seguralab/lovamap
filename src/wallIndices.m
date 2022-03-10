% Find indices of 'edge' voxels along rectangular domain

function voxel_indices = wallIndices(voxelCenters_lin, dx)
    % Store voxels that lie within domain
    % Start by getting x,y,z axes ticks
    zTicks = voxelCenters_lin(1, 3) : dx : voxelCenters_lin(end, 3);
    yTicks = voxelCenters_lin(1, 2) : dx : voxelCenters_lin(end, 2);
    xTicks = voxelCenters_lin(1, 1) : dx : voxelCenters_lin(end, 1);

    % Number of elements in each axis
    shape = [length(xTicks), length(yTicks), length(zTicks)];  

    % Ticks of interest
    edges_log = false(size(yTicks));
    edges_log([1, end]) = true;

    % Get row indices of voxelCenters_lin that correspond to x,y,z
    % For z column, each repeated number is shape(1)*shape(2) long, and
    % there are shape(3) of these repeated numbers
    z_indices = [(1 : shape(1)*shape(2))';
                 ((length(voxelCenters_lin) - shape(1)*shape(2) + 1) : length(voxelCenters_lin))'];
    % For y column, each repeated number is shape(1) long, and there are
    % shape(2) different repeated numbers in a set. There are shape(3)
    % sets.

    y_edges_log = (repelem(edges_log, shape(1)))';
    y_edges_log = repmat(y_edges_log, shape(3), 1);
	y_indices     = (1 : 1 : length(y_edges_log))';
    y_indices     = y_indices(y_edges_log);

    % For x column, there are shape(1) distinct numbers that repeat
    % shape(2) * shape(3) times
    x_edges_log = repmat(edges_log', shape(2) * shape(3), 1);
    x_indices     = (1 : 1 : length(x_edges_log))';
    x_indices     = x_indices(x_edges_log);
    
    % Obtain final indices
%     voxel_indices = fastIntersect(y_indices, x_indices, 'elements');
%     voxel_indices = fastIntersect(voxel_indices, z_indices, 'elements');
    voxel_indices = unique([y_indices; x_indices; z_indices]);

end