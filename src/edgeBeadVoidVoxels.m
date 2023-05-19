% Find the edge voxels of the beads touching edge of void space
 
function edge_beadind = edgeBeadVoidVoxels(in_log, edgeIndices, shape, dx, ...
                                           thickness_scale)
    bin_image = ones(shape);
    bin_image(in_log) = 0;
    % Get EDT from interior to exterior voxels
    if sum(bin_image(:)) == 0 % void space runs to the edge of the image (likely a cropped domain)
        bin_image(1,:,:) = 1; % front
        bin_image(:,1,:) = 1; % left
        bin_image(:,:,1) = 1; % bottom
        bin_image(shape(1),:,:) = 1; % back
        bin_image(:,shape(1),:) = 1; % right
        bin_image(:,:,shape(1)) = 1; % top
        edt_edge3 = find(bin_image);
    else
        edt_edge = dx * bwdist(bin_image);
        edt_edge2 = find(edt_edge <= (dx * thickness_scale * sqrt(2)) + 1e-5);
        edt_edge3 = fastIntersect(edt_edge2, find(in_log), 'elements');
    end
    % store edge voxels per bead
    edge_beadind = cell(length(edgeIndices), 1);
    for i = 1 : length(edgeIndices)
        edge_beadind{i} = fastIntersect(edt_edge3, edgeIndices{i}, 'elements');
    end
end
