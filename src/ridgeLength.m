% Determine arc length of ridges1D by fitting curve
%
% ridge_pts:    indices of ridge points
% ridge_endpts: 1 or 2 peak indices that flank ridge and will be used as anchor
%               points for curve fitting

function lngth = ridgeLength(ridge_pts, ridge_endpts, voxels, dx)

    num_endpts = length(ridge_endpts);

    switch num_endpts
        % single endpoint
        case 1
            % parameterize each dimension
            t = linspace(0, 1, length(ridge_pts) + 1);

            % Peter's sorting function
            ridge_pts_ordered = sortRidgePts(ridge_pts, ridge_endpts, voxels, dx);

            xyz = voxels(ridge_pts_ordered, :);

            [slm_x, ~, x1] = slmengine(t, xyz(:,1), 'leftvalue', voxels(ridge_endpts, 1));
            [slm_y, ~, y1] = slmengine(t, xyz(:,2), 'leftvalue', voxels(ridge_endpts, 2));
            [slm_z, ~, z1] = slmengine(t, xyz(:,3), 'leftvalue', voxels(ridge_endpts, 3));                                     

            lngth = arclength(x1, y1, z1);

        case 2
            % parameterize each dimension
            t = linspace(0, 1, length(ridge_pts) + 2);

            % Peter's sorting function
            ridge_pts_ordered = sortRidgePts(ridge_pts, ridge_endpts, voxels, dx);

            xyz = voxels(ridge_pts_ordered, :);

            [slm_x, ~, x1] = slmengine(t, xyz(:,1), 'leftvalue', voxels(ridge_endpts(1), 1), ...
                                                    'rightvalue', voxels(ridge_endpts(2), 1));
            [slm_y, ~, y1] = slmengine(t, xyz(:,2), 'leftvalue', voxels(ridge_endpts(1), 2), ...
                                                    'rightvalue', voxels(ridge_endpts(2), 2));
            [slm_z, ~, z1] = slmengine(t, xyz(:,3), 'leftvalue', voxels(ridge_endpts(1), 3), ...
                                                    'rightvalue', voxels(ridge_endpts(2), 3));                                     

            lngth = arclength(x1, y1, z1);
        otherwise
            error('Too many or not enough endpoints provided for a ridge.')
    end

%     % Visually check curve fit with 2D plots and 3D plot
%     plotslm(slm_x)
%     plotslm(slm_y)
%     plotslm(slm_z)
%     figure
%     plot3(x1, y1, z1)
%     axis equal
end