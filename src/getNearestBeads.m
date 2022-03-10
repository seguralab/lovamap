%GETNEARESTBEADS LOVAMAP helper function for extracting the indices of the nearest beads to 
% particular voxels from the spares data structure that is currently being used to store information
% of which beads are equidistant to each voxel.
function result = getNearestBeads(sparseData, colIdxs)
    if numel(colIdxs) == 0
        result = [];
    elseif numel(colIdxs) == 1
        [tmp, ~] = find(sparseData(:, colIdxs));
        result = tmp(:);
    else
        % Find the number of nearest beads for each queried voxel
        nNearBeads = full(sum(sparseData(:, colIdxs), 1));
        assert(all(nNearBeads == nNearBeads(1)), ...
               'LOVAMAP:getNearestBeadsError', ...
               'Not all of the queried voxels have the same number of nearest beads.');
        [tmp, ~] = find(sparseData(:, colIdxs));
        result = reshape(tmp, nNearBeads(1), []);
    end
end