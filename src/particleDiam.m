% Return diameter information for particles in domain
% Set chopped-off particles to average particle diameter

function diameters = particleDiam(bead_cell, voxels, dx)
    
    max_z = max(voxels(:, 3));
    diameters = zeros(length(bead_cell), 1);
    adjust = false(length(bead_cell), 1);
    for i = 1 : length(bead_cell)
        % skip if bead intersects with z-axis max, indicating it may
        % not be a full bead
        if binarySearch(voxels(bead_cell{i}, 3), max_z) < 0
            % get volume
            vol = round2sigdig(length(bead_cell{i}) * dx^3, dx);
            % convert to diameter
            diameters(i) = (vol * (3/4) / pi)^(1/3) * 2;
        else
            adjust(i) = true;
        end
    end
    % remove empty diameters
    avg_diam = sum(diameters) / sum(~adjust);
    diameters(adjust) = avg_diam;
end