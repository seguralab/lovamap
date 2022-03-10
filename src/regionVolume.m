function [volume_idxs] = regionVolume(region_idxs, edt, shape, dx)

    % Purposely using a row-vector store my region indices.
    regionIdxs = int32(region_idxs(:)');
    numRegionVoxels = int32(numel(regionIdxs));

    % Casting shape values into int32.
    myShape = int32(shape(:));
    % Strides for the full domain.
    yStride = myShape(1);
    zStride = myShape(1) * myShape(2);

    % regionIJKs is a 3 x numRegionVoxels array where the j-th column stores the
    % ijk of the j-th region voxel index.
    regionIJKs       = zeros(3, numRegionVoxels, 'int32');
    regionIJKs(3, :) = idivide(regionIdxs, zStride) + 1;
    tmpRegionIdxs    = mod(regionIdxs, zStride);
    regionIJKs(2, :) = idivide(tmpRegionIdxs, yStride) + 1;
    regionIJKs(1, :) = mod(tmpRegionIdxs, yStride);

    % Mask for the region volume voxels.
    rvMask = false(int32(prod(shape)), 1);

    % Remember: we are working completely in index space!
    for i = 1 : numRegionVoxels
        ijk = regionIJKs(:, i);
        width = int32(ceil(edt(regionIdxs(i)) / dx));

        % Bounding box for all the voxels we would need to check distances against.
        bboxMin = max(ijk - width, ones(3, 1, 'int32'));
        bboxMax = min(ijk + width, myShape);

        % Cache a list of all the ijks and indices in the bbox.
        % IMPORTANT: Note that mapBBox2World acts as a map from bbox voxel indices to world
        % indices, i.e., if mapBBox2World(i) = j, then the i-th voxel index in the bbox has
        % world index j.
        [ii, jj, kk] = ndgrid(bboxMin(1) : bboxMax(1), ...
                              bboxMin(2) : bboxMax(2), ...
                              bboxMin(3) : bboxMax(3));
        bboxIJKs = [ii(:) jj(:) kk(:)];
        mapBBox2World = bboxIJKs(:, 1) + ...
                       (bboxIJKs(:, 2) - 1) * yStride + ...
                       (bboxIJKs(:, 3) - 1) * zStride;

        % Compute distances.
        d = sum((ijk - bboxIJKs').^2, 1);

        % Identify which distances are small enough to be in the sphere, and map the indices
        % back to world indices so we can record the voxels in the region volume mask.
        rvMask(mapBBox2World(d < width^2)) = true;
    end

    volume_idxs = find(rvMask);
end
