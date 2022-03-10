% Plotting 3D surface
% First input is string indicating whether to show grid or not, i.e., 'on'
%                or 'off'
% Second input:
% 3 possible input formats:
    % Assume varargin = {a, b, c}
    %
    % A) a = 2-column matrix of pixel centers
    %    b = matrix of z-values
    %    c = 'on' or 'off' for showing grid
    %
    % B) a = 3-column matrix of voxel centers
    %    b = matrix of z-values
    %    c = 'on' or 'off' for showing grid
    %
    % C) a = 2-D bounds array: [xmin, xmax, ymin, ymax]
    %    b = dx
    %    c = matrix of z-values
    %    d = 'on' or 'off' for showing grid
    %
    % D) a = 3-D bounds array: [xmin, xmax, ymin, ymax, zmin, zmax]
    %    b = dx
    %    c = matrix of z-values
    %    d = 'on' or 'off' for showing grid
    %

function [] = plotting3Dsurf(varargin)
    
    if nargin < 3
        error('Too few inputs.')
    end
    
    % (A)
    if nargin == 3 && (size(varargin{1}, 2) == 2)
        % a = 2-column matrix of pixel centers
        % b = matrix of z-values
        % c = 'on' or 'off' for showing grid
        
%         if size(varargin{1}, 2) ~= 2
%             error('Pixel centers must be 2-column matrix')
%         end
        if ~ismatrix(varargin{2})
            error('z-values must be a matrix')
        end
        
        % Show grid or not
        if strcmp(varargin{3}, 'off') || strcmp(varargin(3), 'no')
            grid_show = 'none';
        else
            grid_show = 'k';
        end
        
        % *** consider surfl function ***
        surf(reshape(varargin{1}(:,1), size(varargin{2})), ...
             reshape(varargin{1}(:,2), size(varargin{2})), ...
             varargin{2}, 'EdgeColor', grid_show);
    end
    
    % (B)
    if nargin == 3 && (size(varargin{1}, 2) == 3)
        % a = 3-column matrix of voxel centers
        % b = matrix of z-values
        % c = 'on' or 'off' for showing grid
        
        if ~ismatrix(varargin{2})
            error('z-values must be a matrix')
        end
        
        % Show grid or not
        if strcmp(varargin{3}, 'off') || strcmp(varargin(3), 'no')
            grid_show = 'none';
        else
            grid_show = 'k';
        end
        
        % *** consider surfl function ***
        sampleZ = varargin{1}(1,3);
        pixelCenter_lin = varargin{1}(varargin{1}(:,3) == sampleZ, :);
        surf(reshape(pixelCenter_lin(:,1), size(varargin{2})), ...
             reshape(pixelCenter_lin(:,2), size(varargin{2})), ...
             varargin{2}, 'EdgeColor', grid_show);
    end
    
    % (C)
    if nargin == 4 && (size(varargin{1}, 1) == 4)
        % a = bounds array: [xmin, xmax, ymin, ymax]
        % b = dx
        % c = matrix of z-values
        % d = 'on' or 'off' for showing grid
        
        if ~isvector(varargin{1})
            error('Bounds must be a vector')
        end
        if ~isscalar(varargin{2})
            error('dx must be a scalar')
        end
        if ~ismatrix(varargin{3})
            error('z-values must be a matrix')
        end
        
        % Show grid or not
        if strcmp(varargin{3}, 'off') || strcmp(varargin(3), 'no')
            grid_show = 'none';
        else
            grid_show = 'k';
        end
        
        pixelCenter_lin = centeredGrid2D(varargin{1}, varargin{2});
        surf(reshape(pixelCenter_lin(:,1), size(varargin{3})), ...
             reshape(pixelCenter_lin(:,2), size(varargin{3})), ...
             varargin{3}, 'EdgeColor', grid_show);
    end
    
    % (D)
    if nargin == 4 && (size(varargin{1}, 1) == 6)
        % a = bounds array: [xmin, xmax, ymin, ymax, zmin, zmax]
        % b = dx
        % c = matrix of z-values
        % d = 'on' or 'off' for showing grid
        
        if ~isvector(varargin{1})
            error('Bounds must be a vector')
        end
        if ~isscalar(varargin{2})
            error('dx must be a scalar')
        end
        if ~ismatrix(varargin{3})
            error('z-values must be a matrix')
        end
        
        % Show grid or not
        if strcmp(varargin{3}, 'off') || strcmp(varargin(3), 'no')
            grid_show = 'none';
        else
            grid_show = 'k';
        end
        
        voxelCenter_lin = centeredGrid3D(varargin{1}, varargin{2});
        sampleZ = voxelCenter_lin(1,3);
        pixelCenter_lin = voxelCenter_lin(voxelCenter_lin(:,3) == sampleZ, :);
        surf(reshape(pixelCenter_lin(:,1), size(varargin{3})), ...
             reshape(pixelCenter_lin(:,2), size(varargin{3})), ...
             varargin{3}, 'EdgeColor', grid_show);
    end
    
    if nargin > 3
        error('Too many inputs.')
    end

end