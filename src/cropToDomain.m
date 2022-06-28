%CROPTODOMAIN LOVAMAP helper function for cropping beads to a newly cropped domain
function [result] = cropToDomain(eachObj, nVoxels, mask)
    nObjs = numel(eachObj);
    allObjs = int32(cat(1, eachObj{:}));
    gridIds = sparse(allObjs, ones(numel(allObjs), 1, 'int32'), ...
                     repelem(1 : nObjs, cellfun(@numel, eachObj)).', ...
                     nVoxels, 1);
    [ii, ~, ids] = find(gridIds(mask));
    result = accumarray(ids, ii, [], @(r) {sort(r)});
end
