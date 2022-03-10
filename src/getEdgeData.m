% Obtain edge voxels of subunit

% varargin  - (optional) vector of subunit 'center' indices used to compute
%             edge EDT
% varargout - EDT value along edge indices

function [edge_indices, varargout] = getEdgeData(sub_indices, voxels, shape, dx, varargin)

    % find bounding box for each subunit
    domMin = voxels(1, :);
    % shift domain voxels to (1, 1, 1);
    % divide by dx to get matrix voxel as opposed to real-space point
    sub_voxs = voxels(sub_indices, :);
    sub_vox0 = floor((sub_voxs - domMin) / dx) + 1;
    % add 1 voxel buffer
    xmin_i = min(sub_vox0(:, 1)) - 1;
    ymin_i = min(sub_vox0(:, 2)) - 1;
    zmin_i = min(sub_vox0(:, 3)) - 1;
    xmax_i = max(sub_vox0(:, 1)) + 1;
    ymax_i = max(sub_vox0(:, 2)) + 1;
    zmax_i = max(sub_vox0(:, 3)) + 1;

    xIdx   = xmin_i : xmax_i;
    yIdx   = ymin_i : ymax_i;
    zIdx   = zmin_i : zmax_i;
    shape_box = [length(xIdx), length(yIdx), length(zIdx)];

    [xxIdx, yyIdx, zzIdx] = ndgrid(xIdx, yIdx, zIdx);
    allIdx = [xxIdx(:) yyIdx(:) zzIdx(:)];

    allLinIdx   = all_sub2ind(shape, allIdx);
    sub_vox0_bb = sub_vox0 - allIdx(1, :) + 1; % bounding box
    sub_boxind = all_sub2ind(shape_box, sub_vox0_bb);

    % Store edge voxels using bwdist
    bin_imsub = ones(shape_box);
    bin_imsub(sub_boxind) = 0;
    edt_edgesub = dx * bwdist(bin_imsub);
    edt_edge_indsub = find(edt_edgesub <= (dx * sqrt(2)) + 1e-5);
    edge_indsub = intersect(edt_edge_indsub, sub_boxind);
    edge_indices = allLinIdx(edge_indsub);
    
    % if optional input into function is subunit 'center' data
    if ~isempty(varargin) && isvector(varargin{1})
        % Compute EDT from center along edge voxels
        center_voxs = voxels(varargin{1}, :);
        center_vox0 = floor((center_voxs - domMin) / dx) + 1;

        center_vox0_bb = center_vox0 - allIdx(1, :) + 1;
        center_boxind = all_sub2ind(shape_box, center_vox0_bb);

        bin_imsub = zeros(shape_box);
        bin_imsub(center_boxind) = 1;
        edt_centersub = dx * bwdist(bin_imsub);
        varargout = {edt_centersub(edge_indsub)};
    end
end