%CONSTRUCTRIDGES1D LOVAMAP helper function for determining all unique 1D-ridges.
function [ridges1D, orphanMask] = constructRidges1D(medax, loopsBySize)
    medMax   = max(medax.numBeads);
    numBeads = size(medax.beads, 1);

    ridges1D = struct;
    ridges1D.indices_sized = cell(medMax, 1);
    ridges1D.beads_sized   = cell(medMax, 1);
    ridges1D.beads_indexed = cell(medMax, 1);

    orphanMask = true(size(medax.indices));

    for len = cell2mat(loopsBySize.keys)
        % Cache the loops
        loopsOfLen = sortrowsReverse(loopsBySize(len));

        % `len x N` matrix, where each column contains the set of `len` beads
        % that the `N` medax voxels are close to.
        medaxBds = getNearestBeads(medax.beads, medax.sized{len});
        if isempty(medaxBds)
            continue;
        end

        % We want to match each set of all medax voxels' beads to the loops found
        % by sssr, but let's only do it for the unique sets of beads (less work).
        [medaxBdsUniq, ~, ic] = unique(medaxBds.', 'rows');

        % Match the unique sets of medax voxels' beads with the loops, producing the
        % indices of the matches.
        [~, locb] = ismember(medaxBdsUniq, loopsOfLen, 'rows');

        % Use the third output of `unique` to produce the match for all medax voxels.
        medaxMatchLoopIdx = locb(ic);

        % Any 0's in the match indices indicate a medax voxel that was not matched
        % with a loop.
        validLoopIdx = medaxMatchLoopIdx ~= 0;

        % Orphans are the ones that did not get matched.
        % medax.orphans{len} = medax.indices(medax.sized{len}(~validLoopIdx));

        if any(validLoopIdx)
            % Get the unique loop indices that medax voxels were matched to.
            % Each unique loop that was matched produces a new ridge1D.
            % The third output of `unique` is very important! It ultimately
            % tells us which medax voxels got matched to each loop.
            [c, ~, ic] = unique(medaxMatchLoopIdx(validLoopIdx));

            medaxIndTmp = medax.sized{len}(validLoopIdx);

            % Use `ic`, the third output of `unique`, to help us easily partition
            % the medax voxels (via their linear indices) into separate cell arrays,
            % each representing the voxels of a unique ridge1D.
            ridges1D.indices_sized{len} = accumarray(ic, medax.indices(medaxIndTmp), [], @(x) {x});

            % Grab the associated loops of each ridge1D.
            ridges1D.beads_sized{len}   = loopsOfLen(c, :);

            % TODO Deprecate this. But currently, this cell array contains
            % the unique linear ID number of each set of beads that was
            % matched with a ridge1D.
            ridges1D.beads_indexed{len} = all_sub2ind(repelem(numBeads, len), ridges1D.beads_sized{len});

            orphanMask(medax.sized{len}(validLoopIdx)) = false;
        end
    end

    % Number of ridges1D
    ridges1D.num = sum(cellfun(@(x) numel(x), ridges1D.indices_sized));
end