% Peter function for sorting ridge points

function [sortedIdxs] = sortRidgePts(ridgeIdxs, endptIdxs, voxels, dx)
%SORTRIDGEPTS Order the voxels of a ridge including its endpoints so that we may
% find a best-fit curve for the ridge.
%
%  sortedIdxs = SORTRIDGEPTS(ridgeIdxs, endptIdxs, voxels, dx) where
%    ridgeIdxs = 1D array of the flat indices of a ridge
%    endptIdx  = 1D array of the flat indices of the endpoint(s) of the ridge
%    voxels    = N x 3 array of the 3d coordinates of all voxels
%    dx        = voxel size

    % Grab the coordinates of the ridge and endpoints.
    ridgeIdxs(ridgeIdxs < 1) = 1;
    endptIdx(endptIdxs < 1) = 1;
    xyzRidge = voxels(ridgeIdxs, :);
    xyzEnds = voxels(endptIdxs, :);

    % Attach one of the endpoints to the ridge points to create a list of
    % all the nodes.
    nodes = ceil([xyzEnds(1, :); xyzRidge] / dx);

    % Precompute the direction of the 26 neighbors
    nbrs = [...
        -1, -1, -1; ...
        -1, -1,  0; ...
        -1, -1,  1; ...
        -1,  0, -1; ...
        -1,  0,  0; ...
        -1,  0,  1; ...
        -1,  1, -1; ...
        -1,  1,  0; ...
        -1,  1,  1; ...
         0, -1, -1; ...
         0, -1,  0; ...
         0, -1,  1; ...
         0,  0, -1; ...
         % 0,  0,  0; ...
         0,  0,  1; ...
         0,  1, -1; ...
         0,  1,  0; ...
         0,  1,  1; ...
         1, -1, -1; ...
         1, -1,  0; ...
         1, -1,  1; ...
         1,  0, -1; ...
         1,  0,  0; ...
         1,  0,  1; ...
         1,  1, -1; ...
         1,  1,  0; ...
         1,  1,  1; ...
    ];

    % Used to store graph connectivity information
    n1 = [];
    n2 = [];

    for idx = 1 : size(nodes, 1)
        curr = nodes(idx, :);
        tmp = repmat(curr, size(nbrs, 1), 1) + nbrs;
        [~, IA, ~] = intersect(nodes(idx+1 : end, :), tmp, 'rows');
        if ~isempty(IA)
            % Note that because of how we're invoking intersect
            % above, IA is guaranteed to be greater than idx,
            % so we're only saving upper triangular
            % connections.
            n1 = [n1, repmat(idx, 1, numel(IA))];
            n2 = [n2, idx + IA(:)'];
        end

    end

    % List of node indices that are connected, aka, that have neighbors.
    % If there are nodes without neighbors, then those indices will not be
    % in uniqueIdxs.
    uniqueNodeIdxs = unique([n1, n2]);

    % Get the indices of the missing nodes
    missingNodeIdxs = setdiff(1 : size(nodes, 1), uniqueNodeIdxs, 'sorted');

    nodeIdxsInGraph = uniqueNodeIdxs;
    for nodeIdx = missingNodeIdxs
        nns = knnsearch(nodes(nodeIdxsInGraph, :), ...
                        nodes(nodeIdx, :), ...
                        'NSMethod', 'kdtree', ...
                        'IncludeTies', true);
        n1 = [n1, repmat(nodeIdx, 1, length(nns{1}))];
        n2 = [n2, nodeIdxsInGraph(nns{1})];
        nodeIdxsInGraph = [nodeIdxsInGraph, nodeIdx];
    end

    % Construct the graph object
    g = graph(n1, n2);
    % Finally, compute the graph Laplacian - this is what we need! I'm throwing
    % in a negative to keep it consistent with the standard Laplacian operator.
    L = -laplacian(g);

    % Now that we have the graph Laplacian, let's solve the heat equation
    % (diffusion equation) on the graph.
    sz = size(nodes, 1);      % problem size
    h  = 1;                   % time-step for heat/diffusion evolution
    I  = sparse(eye(sz));     % sparse identity matrix

    % Initial condition will be a large amount of heat at the first endpoint.
    uInit = zeros(sz, 1);
    uInit(1) = 100;

    % Propagate the initial heat distribution forward in time by one timestep.
    uSol = (I - h * L) \ uInit;

    % The heat distribution of uSol should be such that the heat at each node
    % decreases as we move farther away from the first endpoint - the endpoint
    % that started with high heat. Obtaining the permutation for sorting this
    % distribution in descending order should tell us the order in which we
    % should order the ridge points.
    [~, permIdxs] = sort(uSol, 'descend');

    % Apply sorted index vector to ridge pts + endpoints for output.
    sortedIdxs = [endptIdxs(1); ridgeIdxs];
    sortedIdxs = sortedIdxs(permIdxs);

    if length(endptIdxs) == 2
        sortedIdxs = [sortedIdxs; endptIdxs(2)];
    elseif length(endptIdxs) > 2
        error('Too many endpoints provided for a ridge.')
    end
end
