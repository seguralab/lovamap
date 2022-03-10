%FINDNEARESTBEADS Summary of this function goes here.
%    Detailed explanation goes here

function [nearestBeads, sparseEDTs, domainEDT] = ...
    findNearestBeads(shape, dx, beadsIdxList, voxelList, tol)

    arguments
        shape (1, 3) int32
        dx (1, 1) double
        beadsIdxList {mustBeCellArray}
        voxelList int32 {mustBeVector}
        tol (1, 1) double = sqrt(2)
    end

    nVoxels = numel(voxelList);
    nBeads = numel(beadsIdxList);
    nDomainVoxels = prod(shape);

    % Some cells of beadsIdxList might be empty, which is fine, but let's
    % construct a vector of indices for the non-empty ondes to use later.
    beadIDs = find(cellfun(@(x) ~isempty(x), beadsIdxList));
    beadIDs = beadIDs(:)';

    domainBBox = CoordBBox(1, shape);

    domainEDT = false(shape);
    domainEDT(cell2mat(beadsIdxList)) = true;
    domainEDT = bwdist(domainEDT); % This is not scaled
    maxValueIS = int32(ceil(max(domainEDT(voxelList))));

    % If worldIdx2VoidIdx(i) = j, then the i-th voxel in the domain is the j-th void voxel.
    worldIdx2VoidIdx = zeros(prod(shape), 1, 'int32');
    worldIdx2VoidIdx(voxelList) = 1 : numel(voxelList);

    % Construct bounding boxes around all beads, simply skipping over the empty beads.
    beadBBoxes = cell(nBeads, 1);
    for i = beadIDs
        [ii, jj, kk] = ind2sub(shape, beadsIdxList{i});
        beadBBoxes{i} = CoordBBox(min([ii(:) jj(:) kk(:)], [], 1), ...
                                  max([ii(:) jj(:) kk(:)], [], 1));
        beadBBoxes{i} = beadBBoxes{i}.expand(maxValueIS).intersect(domainBBox);
    end

    vsMask = false(prod(shape), 1);
    vsMask(voxelList) = true;

    % Must first initialize these cell arrays with empty objects for the case
    % where some of the bead idx lists are empty.
    sparseData = cell(nBeads, 1);
    sparseEDTs = cell(nBeads, 1);
    for i = 1 : nBeads
        sparseData{i} = int32.empty;
        sparseEDTs{i} = sparse([]);
    end

    for i = beadIDs
        [ii, jj, kk] = ind2sub(shape, beadsIdxList{i});
        localBeadIJKs = int32([ii(:) jj(:) kk(:)]) - beadBBoxes{i}.Min + 1;
        localBeadIdxs = all_sub2ind(beadBBoxes{i}.Dim, localBeadIJKs);

        bboxIdxs = all_sub2ind(shape, beadBBoxes{i}.IJKs);

        beadEDT = false(beadBBoxes{i}.Dim);
        beadEDT(localBeadIdxs) = true;
        beadEDT = bwdist(beadEDT);

        nearby = (abs(beadEDT(:) - domainEDT(bboxIdxs))) <= tol;
        nearby = nearby & vsMask(bboxIdxs);

        sparseEDTs{i} = sparse(bboxIdxs(nearby), ...
                               ones(nnz(nearby), 1, 'int32'), ...
                               double(dx * beadEDT(nearby)), ...
                               nDomainVoxels, ...
                               1);

        sparseData{i} = worldIdx2VoidIdx(bboxIdxs(nearby));

    end

    nearestBeads = sparse(repelem(int32(1 : nBeads), cellfun(@numel, sparseData)), ...
                          cell2mat(sparseData), ...
                          true(sum(cellfun(@numel, sparseData)), 1), ...
                          nBeads, ...
                          nVoxels);

    domainEDT = double(dx * domainEDT);
end

% Custon validation function for validating cell arrays
function mustBeCellArray(x)
    if ~iscell(x)
        eid = "findNearestBeads:mustBeCellArray";
        msg = ['The beads index list must be a cell array that contains vectors' ...
               'of the linear indices of voxels that are contained within each bead.'];
        error(msg);
    end
end
