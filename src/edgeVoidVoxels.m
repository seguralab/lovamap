% Find the edge voxels of the void space

function edge_ind = edgeVoidVoxels(bead_voxels, vsVoxels, shape, dx, thickness_scale)
    bin_image = ones(shape);
    bin_image(vsVoxels) = 0;
    bin_image(bead_voxels) = 0;
    % Get EDT from interior to exterior voxels
    edt_edge = dx * bwdist(bin_image);
    edt_edge_ind = find(edt_edge <= (dx * thickness_scale * sqrt(2)) + 1e-5);
    edge_ind = fastIntersect(edt_edge_ind, vsVoxels, 'elements');
end