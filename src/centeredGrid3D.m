% Create 3D grid of centered voxel coordinates

function voxelCenter = centeredGrid3D(bounds, dx)

    xMin = bounds(1);
    xMax = bounds(2);
    yMin = bounds(3);
    yMax = bounds(4);
    zMin = bounds(5);
    zMax = bounds(6);
    
    [xxx, yyy, zzz] = ndgrid(xMin : dx : xMax, ...
                             yMin : dx : yMax, ...
                             zMin : dx : zMax);
                      
    xxx = xxx(1 : end-1, 1 : end-1, 1 : end-1) + (dx / 2);
    yyy = yyy(1 : end-1, 1 : end-1, 1 : end-1) + (dx / 2);
    zzz = zzz(1 : end-1, 1 : end-1, 1 : end-1) + (dx / 2);
    voxelCenter = [xxx(:), yyy(:), zzz(:)];

end