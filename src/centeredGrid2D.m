% Create 2D grid of centered pixel coordinates

function pixelCenter = centeredGrid2D(bounds, dx)

    xMin = bounds(1);
    xMax = bounds(2);
    yMin = bounds(3);
    yMax = bounds(4);
    
    [xx, yy] = meshgrid(xMin : dx : xMax, yMin : dx : yMax);
                      
    xx = xx(1 : end-1, 1 : end-1) + (dx / 2);
    yy = yy(1 : end-1, 1 : end-1) + (dx / 2);
    
    pixelCenter = [xx(:), yy(:)];
end