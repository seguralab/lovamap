% Find the edge voxels of the void space
 
function edge_ind = edgeVoidVoxels(bead_voxels, vsVoxels, shape, dx, thickness_scale)
    bin_image = ones(shape);
    bin_image(vsVoxels) = 0;
    bin_image(bead_voxels) = 0;
    % Get EDT from interior to exterior voxels
    if sum(bin_image(:)) == 0 % void space runs to the edge of the image (likely a cropped domain)
        bin_image(1,:,:) = 1; % front
        bin_image(:,1,:) = 1; % left
        bin_image(:,:,1) = 1; % bottom
        bin_image(shape(1),:,:) = 1; % back
        bin_image(:,shape(1),:) = 1; % right
        bin_image(:,:,shape(1)) = 1; % top
        edge_ind = find(bin_image);
    else
        edt_edge = dx * bwdist(bin_image);
        edt_edge_ind = find(edt_edge <= (dx * thickness_scale * sqrt(2)) + 1e-5);
        edge_ind = fastIntersect(edt_edge_ind, vsVoxels, 'elements');
    end
end
