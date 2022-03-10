% Hildebrand local thickness

function thickness_mat = localThickness(sub_indices, skel_indices, skel_EDT, ...
                                        voxels, dx, shape)
    thickness_mat = zeros(length(sub_indices), 1);
    
    % Domain minimum
    dMin = voxels(1, :);
    
    for i = 1 : length(skel_indices)
        % Locate sub points within bounding box of skeleton point
        % Divide by dx to get matrix voxel as opposed to real-space point
        centerIdx = floor((voxels(skel_indices(i), :) - dMin) / dx) + 1;
        buffer    = ceil(skel_EDT(i) / dx) + 1;
        bboxMin   = max(centerIdx - buffer, [1, 1, 1]);
        bboxMax   = min(centerIdx + buffer, shape);
        xIdx      = bboxMin(1) : bboxMax(1);
        yIdx      = bboxMin(2) : bboxMax(2);
        zIdx      = bboxMin(3) : bboxMax(3);

        [xxIdx, yyIdx, zzIdx] = ndgrid(xIdx, yIdx, zIdx);
        allIdx = [xxIdx(:) yyIdx(:) zzIdx(:)];

        allLinIdx = all_sub2ind(shape, allIdx);
        inds_key = fastIntersect(sub_indices, allLinIdx, 'indices');
        % Distance formula from subunit voxels to current skeleton pt
        d = sum((voxels(skel_indices(i), :) - voxels(sub_indices(inds_key), :)).^2, 2);
        % Logical column vector indicating which sub voxels lie within EDT
        % sphere of current skeleton pt
        these_inds = inds_key(d <= skel_EDT(i).^2);
        
        thickness_mat(these_inds) = max([thickness_mat(these_inds), ...
                                         repmat(skel_EDT(i), numel(these_inds), 1)], ...
                                        [], 2);
    end
    % assign EDT to skeleton indices
    skel_ind_sub = arrayfun(@(x) binarySearch(sub_indices, x), skel_indices);
    thickness_mat(skel_ind_sub) = skel_EDT;
    thickness_mat = thickness_mat * 2;
end