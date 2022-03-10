% medax_ind:          column vector of peak indices
% medax_beads:        [EDT_CACHE_SIZE x num_peaks] matrix
% medax_numBeads:     column vector of number of beads associated with each
%                     peak
% ridges2D_ind:       column vector of ridges2D indices
% ridges2D_beads:     [2 x number of ridges2D]

% medax.indices = Column vector of the linear indices relative to the full domain
%                 of the voxels of the medial axis.
%
% ridges2D.cum_indices = Column vector of the linear indices relative to the full
%                        domain of the voxels of 2D ridges.

function [medax, ridges2D] = determineVoxelType(nearestBeads, vsVoxelIdxs)
    medax = struct;
    ridges2D = struct;

    numNearBeads = full(sum(nearestBeads, 1));

    medaxMask = numNearBeads >= 3;
    medax.indices = vsVoxelIdxs(medaxMask);
    medax.beads = nearestBeads(:, medaxMask);
    medax.numBeads = numNearBeads(medaxMask)'; % The transpose is necessary

    ridges2DMask = numNearBeads == 2;
    ridges2D.cum_indices = vsVoxelIdxs(ridges2DMask);
    [beadIdx, ~] = find(nearestBeads(:, ridges2DMask));
    ridges2D.cum_beads = reshape(beadIdx, 2, numel(ridges2D.cum_indices));
end
