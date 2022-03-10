%FINDDISTTOCLOSESTFARBEAD LOVAMAP helper function for finding the distance to the nearest
% bead given a list of beads to ignore.
function [dist, closestVoxelIdxs] = findDistToClosestFarBead(data, beadsToIgnore, voxelList)
    % Create array of indices of beads that does not include the beads that 
    % are close to the ridge1D
    beadIdxs = 1 : numel(data.beads);
    beadIdxs(beadsToIgnore) = []; % Remove indices of beads that are close

    % Create EDT
    localEdt = false(data.shape);
    localEdt(cat(1, data.beads{beadIdxs})) = true;
    [localEdt, closestIdx] = bwdist(localEdt); % Note that this is unscaled by dx.

    dist = data.dx * localEdt(voxelList) - data.EDT(voxelList);
    closestVoxelIdxs = closestIdx(voxelList);
end