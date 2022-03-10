function [domain, dx, res] = setVoxelSize(domain, dxInput, voxelRange)
 
    if numel(voxelRange) == 1
        voxelRange(2) = voxelRange(1);
    end
 
    xExtent = domain(2) - domain(1);
    yExtent = domain(4) - domain(3);
    zExtent = domain(6) - domain(5);
    totalVoxels = xExtent * yExtent * zExtent;
 
    dxMin = (totalVoxels / voxelRange(2)) ^ (1/3);
    dxMax = (totalVoxels / voxelRange(1)) ^ (1/3);
 
    dxIdeal = min(max(dxInput, dxMin), dxMax);
    dx = round2sigdig_exact(dxIdeal, dxInput, 'down');
    if dx == 0
        dx = round2sigdig_exact(dxIdeal, 0.5, 'down');
    end
 
    nx = ceil(xExtent / dx);
    ny = ceil(yExtent / dx);
    nz = ceil(zExtent / dx);
 
    domain(2) = domain(1) + nx * dx;
    domain(4) = domain(3) + ny * dx;
    domain(6) = domain(5) + nz * dx;
 
    res = nx * ny * nz;
end