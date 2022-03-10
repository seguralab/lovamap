% Return a bead struct containing the following fields:
% Beads       : cell array containing indices of each bead
% EdgeIndices : cell array containing edge indices of each bead
% Shell       : cell array containing a 'shell' around beads of given
%               thickness
% Shape       : domain shape
% DomainMin   : domain minimum
% AllBeads    : column vector containing indices of all beads
% AllEdges    : column vector containing edge indices of all beads

function bead_struct = labelBeadDomain(bead_data, voxels, shell_thickness, dx, shape)

    arguments
        bead_data       (:, 4) double
        voxels          (:, 3) double
        shell_thickness (1, 1) double {mustBePositive}
        dx              (1, 1) double {mustBePositive}
        shape           (1, 3) int32  {mustBePositive}
    end

    % Domain minimum
    dMin = voxels(1, :);
    % Number of beads
    num_beads = size(bead_data, 1);
    % Convert shell_thickness to number of voxels
    shell_voxs = shell_thickness / dx;

    domainBBox = CoordBBox(1, shape);

    % Initialize bead cell array
    bead_struct.Beads       = cell(num_beads, 1);
    bead_struct.EdgeIndices = cell(num_beads, 1);
    bead_struct.Shell       = cell(num_beads, 1);

    % Mask of the domain where true means no bead has claimed the voxel yet.
    domainMask = true(shape);

    for i = 1 : num_beads
        centerXYZ = bead_data(i, 1 : 3);
        radius    = bead_data(i, 4);
        centerIJK = floor((centerXYZ - dMin)/ dx) + 1;

        width = ceil(bead_data(i, 4) / dx) + 1;
        bbox = CoordBBox(centerIJK).expand(width).intersect(domainBBox);

        if bbox.isEmpty()
            bead_struct.Beads{i} = int32.empty;
            bead_struct.EdgeIndices{i} = int32.empty;
            bead_struct.Shell{i} = int32.empty;
            continue;
        end

        bboxIdxs = all_sub2ind(shape, bbox.IJKs); % Vertical array

        d = sum((centerXYZ - voxels(bboxIdxs, 1:3)).^2, 2); % Vertical array

        % Mask of bead voxels (size of the bounding box)
        beadMask = (d <= radius^2); % Vertical array

        % Intersect bead mask with a mask of the voxels in this bounding box
        % that are still in the domain, i.e., haven't been claimed by another bead.
        beadMask = beadMask & domainMask(bboxIdxs); % Vertical array

        % Update domain mask to by removing the voxels claimed by this bead.
        domainMask(bboxIdxs) = domainMask(bboxIdxs) & ~beadMask;

        % Store bead indices in cell array
        bead_struct.Beads{i} = bboxIdxs(beadMask);

        % Store edge voxels using bwdist
        edtEdge = true(bbox.Dim);
        edtEdge(beadMask) = false;
        edtEdge = bwdist(edtEdge); % Not scaled by dx
        edtEdgeMask = (edtEdge(:) <= sqrt(2));
        bead_struct.EdgeIndices{i} = bboxIdxs(edtEdgeMask & beadMask);

        % Store a shell cloud of points that is shell_thickness thick
        edtShellMask = (edtEdge(:) <= shell_voxs * sqrt(2));
        bead_struct.Shell{i} = bboxIdxs(edtShellMask & beadMask);
    end

    % Additional information for struct
    bead_struct.Shape     = shape;
    bead_struct.DomainMin = dMin;

    % Get cumulative bead indices list
    bead_struct.AllBeads = cat(1, bead_struct.Beads{:});
    bead_struct.AllEdges = cat(1, bead_struct.EdgeIndices{:});
    bead_struct.AllShells = cat(1, bead_struct.Shell{:});
end
