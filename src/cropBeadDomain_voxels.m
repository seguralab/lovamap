% Crop labeled bead domain
% voxels_cropped: updated voxels
% domain_cropped: xmin...zmax of cropped domain
% voxel_log:      indices of rows from original voxels

function [voxels_cropped, domain_cropped, voxel_log] = ...
                     cropBeadDomain_voxels(voxels, domain, cropPercent, dx)

    nVPDx = (domain(2) - domain(1)) / dx;
    nVPDy = (domain(4) - domain(3)) / dx;
    nVPDz = (domain(6) - domain(5)) / dx;

    a = 1 - cropPercent;
    nVoxelsToTrimX = ceil(nVPDx * a);
    nVoxelsToTrimY = ceil(nVPDy * a);
    nVoxelsToTrimZ = ceil(nVPDz * a);

    if mod(nVoxelsToTrimX, 2) == 1
        nVoxelsToTrimX = nVoxelsToTrimX + 1;
    end
    if mod(nVoxelsToTrimY, 2) == 1
        nVoxelsToTrimY = nVoxelsToTrimY + 1;
    end
    if mod(nVoxelsToTrimZ, 2) == 1
        nVoxelsToTrimZ = nVoxelsToTrimZ + 1;
    end

    portionTrim = dx * 0.5 * [nVoxelsToTrimX, nVoxelsToTrimY, nVoxelsToTrimZ];

    xmin = domain(1) + portionTrim(1);
    xmax = domain(2) - portionTrim(1);
    ymin = domain(3) + portionTrim(2);
    ymax = domain(4) - portionTrim(2);
    % Determine domain to plot
    if (domain(6) - domain(5)) <= ...
       ((domain(2) - domain(1)) * cropPercent)

        zmin = domain(5);
        zmax = domain(6);
    else
        zmin = domain(5) + portionTrim(3);
        zmax = domain(6) - portionTrim(3);
    end
    domain_cropped = [xmin, xmax, ymin, ymax, zmin, zmax];

    % Store voxels that lie within domain
    % Start by getting x,y,z axes ticks
    zTicks = voxels(1, 3) : dx : voxels(end, 3);
    yTicks = voxels(1, 2) : dx : voxels(end, 2);
    xTicks = voxels(1, 1) : dx : voxels(end, 1);

    % Number of elements in each axis
    shape = [length(xTicks), length(yTicks), length(zTicks)];

    % Crop along these ticks and store indices
    zcrop = (zTicks >= zmin) & (zTicks <= zmax);
    zcrop_idx = find(zcrop == 1);

    ycrop = (yTicks >= ymin) & (yTicks <= ymax);

    xcrop = (xTicks >= xmin) & (xTicks <= xmax);
    % Get row indices of voxels that correspond to x,y,z crop
    % For z column, each repeated number is shape(1)*shape(2) long, and
    % there are shape(3) of these repeated numbers
    zminZmax_idx = [1 + (zcrop_idx(1) - 1) * (shape(1)*shape(2)); ...
                            zcrop_idx(end) * (shape(1)*shape(2))];
    z_croplogical = false(size(voxels, 1), 1);
    z_croplogical(zminZmax_idx(1) : zminZmax_idx(2)) = true;

    % For y column, each repeated number is shape(1) long, and there are
    % shape(2) different repeated numbers in a set. There are shape(3)
    % sets.
    y_croplogical = (repelem(ycrop, shape(1)))';
    y_croplogical = repmat(y_croplogical, shape(3), 1);

    % For x column, there are shape(1) distinct numbers that repeat
    % shape(2) * shape(3) times
    x_croplogical = repmat(xcrop', shape(2) * shape(3), 1);

    % Obtain final indices
    voxel_log      = x_croplogical & y_croplogical & z_croplogical;
    voxels_cropped = voxels(voxel_log, :);
end
