% Heat map of 'ligand proteins' in domain (both throughout void space
% and within beads). We assume the RGD concentration provided is
% distributed throughout the entire bead.
%
%
% RGD_location   : 'shell'                only consider RGD in the shell of
%                                         the bead (0 everywhere else)
%                  'full' or 'bead'       consider RGD throughout entire bead
% voxel_datatype : 'cell' or 'average'    assumes a cell is moving around the space and
%                                         sensing the AVERAGE ligand concentration stored at each voxel
%                  'molecules' or 'exact' stores the number of ligand molecules at each
%                                         voxel (only bead voxels contain ligand)
% filter_type    : 'cube'                 use a cubice convolution filter
%                  'sphere'               use a spherical convolution filter
% output_type    : 'all'   output ligand information for entire domain
%                : 'void'  output ligand information for void space only
%                          and set bead voxels to 0
%                : 'bead'  output ligand information for bead only and
%                          set void space voxels to 0
% bead_struct    : contains Beads field (cell array of voxel indices of beads)
%                  and Shell field (cell array of voxel indices of bead shell)
% cell_diameter  : in um (or units) - used to determine the length of the
%                  'sliding window' cube used to form the ligand density
%                  map
% shell_thickness: number of voxels thick
% RGD_conc       : RGD concentration in micromoles/liter per particle
% vsVoxel_log    : boolean vector indicating void space voxel or not
% shape          : domain shape, [xDim, yDim, zDim]
%
% ligand_map     : 3D matrix representing ligand density within domain

function ligand_map = ligandMap(RGD_location, voxel_datatype, filter_type, output_type, ...
                                    bead_struct, cell_diameter, ...
                                    RGD_conc, vsVoxel_log, shape, dx)

    if ~ismember(RGD_location, {'shell', 'full', 'bead'})
        error('RGD_location must be ''shell'', ''full'', or ''bead''.');
    end

    if ~ismember(voxel_datatype, {'cell', 'average', 'molecule', 'molecules', 'exact'})
        error('voxel_datatype must be ''cell'', ''average'', ''molecule(s)'', or ''exact''.')
    end

    % Saving time / memory if these are the inputs...
    if strcmp(output_type, 'void')
        if ismember(voxel_datatype, {'molecule', 'molecules', 'exact'})
            ligand_map = zeros(shape);
            fprintf('%30s %s\n', ...
                'Zero ligand protein concentration was outputted because output_type ''void'' and voxel_datatype ''molecules'' together negate all data.');
            return;
        end
    elseif strcmp(output_type, 'bead')
        beadMask = false(prod(shape), 1);
        beadMask(bead_struct.AllBeads) = true;
    end

    % ensure cell diameter is an integer
    cell_diam_voxs = floor(cell_diameter / dx);
    % ensure cell diameter window is an odd number
    if mod(cell_diam_voxs, 2) < 1
        cell_diam_voxs = cell_diam_voxs - 1;
    end

    % Initialize initial ligand concentration
    ligandInit = zeros(shape);

    if ismember(RGD_location, {'bead', 'full'})
        % Compute number of ligand molecules in bead
        % Using conversions:
        % 1 um^3 = 1e-15 L
        % 1 micromole = 6.02e17 molecules
        nMolsPerVoxel = RGD_conc * dx^3 * 602; % 602 = 1e-15 * 6.02e17;

        ligandInit(data.allBeads) = nMolsPerVoxel;
    else
        % TODO - precompute the concatenation of shell indices
        allShell = unique(cat(1, bead_struct.Shell{:}));

%         shellVols = beadVols - ((1 / 6) * pi * ((6 * beadVols / pi) .^ (1/3) - ...
%                     2 * shell_thickness) .^ 3); % in um^3

        nMolsPerVoxel = RGD_conc * dx^3 * 602; % 602 = 1e-15 * 6.02e17;
        ligandInit(allShell) = nMolsPerVoxel;
    end

    % Populate output matrix based on desired data type
    if ismember(voxel_datatype, {'molecule', 'molecules', 'exact'})
        % store number of molecules at each voxel
        ligand_map = ligandInit;
    else % average is default

        % Cube and sphere filter are the only options
        if strcmp(filter_type, 'cube')
            filter = ones(cell_diam_voxs, cell_diam_voxs, cell_diam_voxs);
        else
            filter = zeros(cell_diam_voxs, cell_diam_voxs, cell_diam_voxs);
            c = ceil(cell_diam_voxs / 2);
            r = floor(cell_diam_voxs / 2);
            rng = 1 : cell_diam_voxs;
            [x, y, z] = meshgrid(rng, rng, rng);
            filter((x-c).^2 + (y-c).^2 + (z-c).^2 <= r^2) = 1;
        end

        ligand_map = convn(ligandInit, filter, 'same') ./ ...
                      convn(ones(size(ligandInit)), filter, 'same');
    end

    % Only show ligand values for desired voxels (0 out all other voxels)
    if strcmp(output_type, 'void')
        if ismember(voxel_datatype, {'molecule', 'molecules', 'exact'})
            % We already checked this invalid case above
        else
            ligand_map(~vsVoxel_log) = 0;
        end
    elseif strcmp(output_type, 'bead')
        if ismember(voxel_datatype, {'molecule', 'molecules', 'exact'})
            % Integrin map is correct
        else
            ligand_map(~beadMask) = 0;
        end
    else % Default is all
    	% Integrin map is correct
    end
end
