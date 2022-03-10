%FINDLOOPS LOVAMAP helper function that finds the smallest set of smallest loops from a bead graph.
function [loopsBySize] = findLoops(medax, beadNeighbors)
    globalLoopMin = min(medax.numBeads);
    globalLoopMax = max(medax.numBeads);
    numBeads = size(medax.beads, 1);

    remainingMask = true(1, numel(medax.indices));
    loops = cell(numBeads, 1);

    for i = 1 : numBeads
        % Subgraph will contain all loops with a minimum node index of i.
        % This mask tell us which voxels to use for its "nearest bead" data
        % to construct our subgraph.
        maskForSubgraph = remainingMask & medax.beads(i, :);
        if all(maskForSubgraph == false)
            continue;
        end

        maxLoopLength = max(medax.numBeads(maskForSubgraph));
        [bdIdxs, ~] = find(medax.beads(:, maskForSubgraph));
        bdIdxs = unique(bdIdxs);

        % Extract the rows of ridges2D.beadsMatrix where both components are
        % members of bdIdxs. These rows make up the subgraph that we're
        % about to process with sssr.
        whichRows = all(ismember(beadNeighbors, bdIdxs), 2);
        tmpLoops = sssr(beadNeighbors(whichRows, :), maxLoopLength);

        % Remove loops that do not contain i.
        tmpLoops = tmpLoops(any(tmpLoops == i, 2), :);

        % Early out if there are no more loops left
        if isempty(tmpLoops), continue; end

         % j tells us the length of each loop
        [~, j] = find(tmpLoops(:,2:end) == tmpLoops(:,1));
        u = unique(j);

        % Convert data into a map, mapping loop length to the set of all such loops.
        sepLoops = arrayfun(@(n) tmpLoops(j==n,1:n), u, 'UniformOutput', false);
        if numel(sepLoops) == 1
            loops{i} = containers.Map(u, sepLoops{:});
        else
            loops{i} = containers.Map(u, sepLoops);
        end

        % Update the mask, turning off the voxels that we've already
        % used the "nearest bead" data of to construct our subgraph.
        remainingMask(medax.beads(i, :)) = false;
    end

    % Remove any cells that are empty.
    loops = loops(~cellfun(@isempty, loops));

    % Convert all maps into a single map, mapping loop length to the set of all such loops.
    loopsBySize = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for len = globalLoopMin : globalLoopMax
        hasLen = cellfun(@(x) x.isKey(len), loops);
        if any(hasLen == true, 'all')
            loopNodes = cellfun(@(x) x(len), loops(hasLen), 'UniformOutput', false);
            loopsBySize(len) = sortrows(sort(vertcat(loopNodes{:}), 2));
        end
    end
end