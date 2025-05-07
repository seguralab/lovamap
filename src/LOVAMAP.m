%LOVAMAP This function analyzes the void space of 3D-packed particles.
function [data, time_log] = LOVAMAP(domain_file, voxel_size, voxel_range, crop_percent, ...
    dip_percent, hall_cutoff, shell_thickness, num_2D_slices, combine_edge_subs)

    arguments
        domain_file {mustBeFile}
        voxel_size (1, 1) double {mustBePositive, mustBeFinite, mustBeNonNan}
        voxel_range (1, 2) double {mustBePositive, mustBeFinite, mustBeNonNan}
        crop_percent (1, 1) double {mustBeInRange(crop_percent, 0, 1)}
        dip_percent (1, 1) double {mustBeNonnegative, mustBeFinite, mustBeNonNan}
        hall_cutoff (1, 1) double
        shell_thickness (1, 1) double {mustBePositive, mustBeFinite, mustBeNonNan}
        num_2D_slices (1, 1) int32 {mustBeInteger}
        combine_edge_subs (1, 1) logical = true;
    end

    totalTimeStart = tic;

    dx = voxel_size;
    tot_time = 0;

    % Struct used to store output data
    data = struct;

    % Struct array for logging time spent in each step of LOVAMAP
    time_log = struct();
    timeLogIdx = length(time_log);

    % Save the input arguments into the output struct for easy reference later
    data.Timestamp                 = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    data.InputArgs                 = struct;
    data.InputArgs.DomainFile      = GetFullPath(domain_file);
    data.InputArgs.VoxelSize       = voxel_size;
    data.InputArgs.VoxelRangeMin   = voxel_range(1);
    data.InputArgs.VoxelRangeMax   = voxel_range(2);
    data.InputArgs.CropPercent     = crop_percent;
    data.InputArgs.DipPercent      = dip_percent;
    data.InputArgs.HallCutoff      = hall_cutoff;
    data.InputArgs.ShellThickness  = shell_thickness;
    data.InputArgs.Num2DSlices     = num_2D_slices;
    data.InputArgs.CombineEdgeSubs = combine_edge_subs;

    % Check file type and read in data
    tStart = tic;
    [~, ~, fExt] = fileparts(domain_file);
    switch lower(fExt)
        case '.json'
            % labeled bead domain
            file = fopen(domain_file);
            raw = fread(file,inf);
            str = char(raw');
            fclose(file);
            val = jsondecode(str);
            % get struct field names
            struct_names = fieldnames(val);
            bead_numbers = fieldnames(val.(struct_names{end}));
            % domain range
            for i = 1 : length(struct_names)
                if strcmp(struct_names{i}, 'domain_size')
                    shape = val.(struct_names{i});
                    shape = shape(:)';
                    domain = [0, shape(1), 0, shape(2), 0, shape(3)];
                    break;
                end
            end
            % dx
            for i = 1 : length(struct_names)
                if strcmp(struct_names{i}, 'voxel_size')
                    dx = val.(struct_names{i});
                    break;
                end
            end
            % number of beads
            for i = 1 : length(struct_names)
                if strcmp(struct_names{i}, 'bead_count')
                    num_beads = val.(struct_names{i});
                    if num_beads ~= length(bead_numbers)
                        error('.json data file: bead_count does not match number of beads listed')
                    end
                    break;
                end
            end

            % number of voxels and beads
            % Shape
            nVPDx = (domain(2) - domain(1)) / dx;
            nVPDy = (domain(4) - domain(3)) / dx;
            nVPDz = (domain(6) - domain(5)) / dx;

            if abs(round(nVPDx) - nVPDx) < 1e-12
                nVPDx = round(nVPDx);
            else
                error('Number of voxels in x is not an integer.')
            end
            if abs(round(nVPDy) - nVPDy) < 1e-12
                nVPDy = round(nVPDy);
            else
                error('Number of voxels in y is not an integer.')
            end
            if abs(round(nVPDz) - nVPDz) < 1e-12
                nVPDz = round(nVPDz);
            else
                error('Number of voxels in z is not an integer.')
            end

            shape = [nVPDx, nVPDy, nVPDz];
            nVoxels = nVPDx * nVPDy * nVPDz;

            % Create 3D grid of voxel coordinates (center of voxel cube)
            voxels = centeredGrid3D(domain, dx);

            % get bead data
            bead_struct = struct;
            bead_struct.Beads = cell(num_beads, 1);
            remove_beads = false(num_beads, 1);
            for i = 1 : num_beads
                ranges = val.(struct_names{end}).(bead_numbers{i});
                bead_struct.Beads{i} = cell2mat(arrayfun(@(x) (ranges(x, 1):ranges(x, end))', ...
                                                         1 : size(ranges, 1), 'UniformOutput', false)');
                if isempty(bead_struct.Beads{i})
                    remove_beads(i) = true;
                end
            end
            bead_struct.Beads(remove_beads) = [];
            num_beads = num_beads - sum(remove_beads);

            if crop_percent < 1
                nVoxelsOld = nVoxels;

                [voxels, domain, crop_mask] = cropBeadDomain_voxels(voxels, domain, ...
                                                         crop_percent, dx);
                % Shape
                nVPDx = (domain(2) - domain(1)) / dx;
                nVPDy = (domain(4) - domain(3)) / dx;
                nVPDz = (domain(6) - domain(5)) / dx;
                %shape_cropped = [nVPDx, nVPDy, nVPDz];

                if abs(round(nVPDx) - nVPDx) < 1e-12
                    nVPDx = round(nVPDx);
                else
                    error('Number of voxels in x is not an integer.')
                end
                if abs(round(nVPDy) - nVPDy) < 1e-12
                    nVPDy = round(nVPDy);
                else
                    error('Number of voxels in y is not an integer.')
                end
                if abs(round(nVPDz) - nVPDz) < 1e-12
                    nVPDz = round(nVPDz);
                else
                    error('Number of voxels in z is not an integer.')
                end
                nVoxels = nVPDx * nVPDy * nVPDz;
                shape = [nVPDx, nVPDy, nVPDz];
                % Update beads
                bead_struct.Beads = cropToDomain(bead_struct.Beads, ...
                                                 nVoxelsOld, crop_mask);
                remove_beads = false(size(bead_struct.Beads, 1), 1);
                for i = 1 : size(bead_struct.Beads, 1)
                    if isempty(bead_struct.Beads{i})
                        remove_beads(i) = true;
                    end
                end
                bead_struct.Beads(remove_beads) = [];
                num_beads = size(bead_struct.Beads, 1);
            end

            data.numVoxels = nVoxels;
            data.dx = dx;
            data.domain = domain;
            data.shape = shape;
            % fprintf('%30s %1.1e\n', 'Number of voxels:', nVoxels);
            % fprintf(runtimes_file, '%30s %1.1e\n', 'Number of voxels:', nVoxels);
            % fprintf('%30s %i\n', 'Number of beads:', num_beads);
            % fprintf(runtimes_file, '%30s %i\n', 'Number of beads:', num_beads);

            % Additional information needed to store full bead data
            domMin     = voxels(1, :);
            shell_voxs = shell_thickness / dx;
            bead_struct.EdgeIndices = cell(num_beads, 1);
            bead_struct.Shell       = cell(num_beads, 1);
            bead_struct.AllBeads    = [];
            bead_struct.AllEdges    = [];

            bead_centers = zeros(num_beads, 3);
            for i = 1 : num_beads
                % Find bead centers to compute convex hull
                bead_centers(i, :) = mean(voxels(bead_struct.Beads{i}, :), 1);

                % Additional bead information
                % Store edge voxels using bwdist
                % shift domain voxels to (1, 1, 1);
                % divide by dx to get matrix voxel as opposed to real-space point
                bead_voxs = voxels(bead_struct.Beads{i}, :);
                bead_vox0 = floor((bead_voxs - domMin) / dx) + 1;
                % add 1 voxel buffer
                xmin_i = min(bead_vox0(:, 1)) - 1;
                ymin_i = min(bead_vox0(:, 2)) - 1;
                zmin_i = min(bead_vox0(:, 3)) - 1;
                xmax_i = max(bead_vox0(:, 1)) + 1;
                ymax_i = max(bead_vox0(:, 2)) + 1;
                zmax_i = max(bead_vox0(:, 3)) + 1;

                xIdx   = xmin_i : xmax_i;
                yIdx   = ymin_i : ymax_i;
                zIdx   = zmin_i : zmax_i;
                shape_box = [length(xIdx), length(yIdx), length(zIdx)];

                [xxIdx, yyIdx, zzIdx] = ndgrid(xIdx, yIdx, zIdx);
                allIdx = [xxIdx(:) yyIdx(:) zzIdx(:)];

                allLinIdx         = all_sub2ind(shape, allIdx);
                bead_vox0_bb      = bead_vox0 - allIdx(1, :) + 1; % bounding box
                bead_boxind       = all_sub2ind(shape_box, bead_vox0_bb);

                % Store edge voxels using bwdist
                bin_imbead = ones(shape_box);
                bin_imbead(bead_boxind) = 0;
                edt_edgebead = dx * bwdist(bin_imbead);
                edt_edge_indbead = find(edt_edgebead <= (dx * sqrt(2)) + 1e-5);
                edge_indbead = intersect(edt_edge_indbead, bead_boxind);
                bead_struct.EdgeIndices{i} = allLinIdx(edge_indbead);

                % Store a shell cloud of points that is shell_thickness thick
                edt_shell_ind = find(edt_edgebead <= (shell_voxs * dx * sqrt(2)) + 1e-5);
                shell_ind = intersect(edt_shell_ind, bead_boxind);
                bead_struct.Shell{i} = allLinIdx(shell_ind);
            end

            % Additional information for struct
            bead_struct.Shape     = shape;
            bead_struct.DomainMin = domMin;
            bead_struct.AllEdges = unique(bead_struct.AllEdges);

            % Compute particle diameter distribution
            diameters = particleDiam(bead_struct.Beads, voxels, dx);
            bead_diams = diameters;

            % Store
            data.beads        = bead_struct.Beads;
            data.edgeIndices  = bead_struct.EdgeIndices;
            data.shell        = bead_struct.Shell;
            data.allBeads     = cat(1, data.beads{:});
            data.allEdges     = unique(cat(1, data.edgeIndices{:}));

        case {'.dat', '.txt'}
            % Read file and convert to array variable
            % Skip hashtags
            fid = fopen(domain_file);
            row_mrkr = fgetl(fid);
            row_cntr = 1;
            data_type = 'default';
            while row_mrkr(1) == '#'
                % Check for data type comment
                if strcmp(row_mrkr(1 : 13), '# Data type: ')
                    if strcmp(sscanf(row_mrkr(14 : end), '%s'), 'Spherical')
                        data_type = 'spherical';
                    elseif strcmp(sscanf(row_mrkr(14 : end), '%s'), 'Labeled')
                        data_type = 'labeled';
                    end
                elseif strcmp(row_mrkr(1 : 14), '# Voxel size: ')
                    dx = sscanf(row_mrkr(15 : end), '%f');
                elseif strcmp(row_mrkr(1 : 15), '# Domain size: ')
                    row_mrkr(strfind(row_mrkr, 'x')) = [];
                    domain_maxes = sscanf(row_mrkr(16 : end), '%f')';
                    % create domain variable that includes domain mins
                    domain = [0, domain_maxes(1), 0, domain_maxes(2), 0, domain_maxes(3)];
                    clear('domain_maxes');
                end
                row_mrkr = fgetl(fid);
                row_cntr = row_cntr + 1;
            end
            fclose(fid);
            % For now, default will be old spherical data
            old_format = false;
            if strcmp(data_type, 'default')
                old_format = true;
                data_type  = 'spherical';
            end
            % Determine datatype
            switch data_type
                case 'spherical'
                    % Read in file
                    if old_format
                        sep = ' ';
                    else
                        sep = ',';
                    end
                    beads = dlmread(domain_file, sep, row_cntr, 0);

                    % Columns 1:3 = (x,y,z) of particle center, Column 4 = particle radii
                    bead_data = beads(:, 1:4);

                   %%%%%%%%%%%%%%%%%% !!!!!!!!!!!!!!!!!!!
                   %%% REMOVE BEADS THAT LIE ABOVE z = 600 %%%
                   rmv_beads = (bead_data(:, 3) + bead_data(:, 4)) > 600;
                   bead_data(rmv_beads, :) = [];

                    % Scan the beads to find min/max bounds and max radius
                    a = min(beads(:, 1:3));
                    b = max(beads(:, 1:3));
                    rMax = max(bead_data(:, end));

                    domain = [a(1) - rMax, b(1) + rMax, ...
                              a(2) - rMax, b(2) + rMax, ...
                              a(3) - rMax, b(3) + rMax];

                    if crop_percent < 1
                        % Crop beads by only including cropPercent of domain
                        [bead_data, domain] = cropBeadDomain(bead_data, domain, crop_percent);
                    end

                    % Set resolution of spheres
                    sphere_res = 30;
                    data.beadDiamLargest = max(bead_data(:, 4)) * 2;
                    data.beads_dat = bead_data;

                    %%%%%%%************************************%%%%%%%
                    %%%%%***************** GRID *****************%%%%%
                    %%%%%%%************************************%%%%%%%
                    % Relabeling input for readability
                    % mins
                    for i = [1 3 5]
                        if domain(i) ~= 0
                            domain(i) = round2sigdig(domain(i), dx, 'down');
                        end
                    end

                    if ~isempty(voxel_range)
                        % Petes code to determine dx
                        [domain, dx, ~] = setVoxelSize(domain, dx, voxel_range);
                        if dx == 0
                            % need significant digit of dx to be decimal
                            [domain, dx, ~] = setVoxelSize(domain, 0.5, voxel_range);
                        end
                    end

                    % fprintf('%29s %.1f\n\n', 'dx:', dx);
                    % fprintf(runtimes_file, '%29s %.1f\n\n', 'dx:', dx);

                    % Create 3D grid of voxel coordinates (center of voxel cube)
                    voxels = centeredGrid3D(domain, dx);

                    % Shape
                    nVPDx = (domain(2) - domain(1)) / dx;
                    nVPDy = (domain(4) - domain(3)) / dx;
                    nVPDz = (domain(6) - domain(5)) / dx;

                    if abs(round(nVPDx) - nVPDx) < 1e-12
                        nVPDx = round(nVPDx);
                    else
                        error('Number of voxels in x is not an integer.')
                    end
                    if abs(round(nVPDy) - nVPDy) < 1e-12
                        nVPDy = round(nVPDy);
                    else
                        error('Number of voxels in y is not an integer.')
                    end
                    if abs(round(nVPDz) - nVPDz) < 1e-12
                        nVPDz = round(nVPDz);
                    else
                        error('Number of voxels in z is not an integer.')
                    end
                    nVoxels = nVPDx * nVPDy * nVPDz;

                    shape = [nVPDx, nVPDy, nVPDz];
                    data.numVoxels = nVoxels;
                    data.dx = dx;
                    data.domain = domain;
                    data.shape = shape;

                    % Transform sphere data into bead struct
                    bead_struct = labelBeadDomain(bead_data, voxels, ...
                                                  shell_thickness, dx, shape);

                    if crop_percent < 1
                        % nVoxelsOld = nVoxels;
                        %
                        % [voxels, domain, crop_mask] = cropBeadDomain_voxels(voxels, domain, ...
                        %                                             crop_percent, dx);
                        % % Shape
                        % nVPDx = (domain(2) - domain(1)) / dx;
                        % nVPDy = (domain(4) - domain(3)) / dx;
                        % nVPDz = (domain(6) - domain(5)) / dx;
                        %
                        % if abs(round(nVPDx) - nVPDx) < 1e-12
                        %     nVPDx = round(nVPDx);
                        % else
                        %     error('Number of voxels in x is not an integer.')
                        % end
                        % if abs(round(nVPDy) - nVPDy) < 1e-12
                        %     nVPDy = round(nVPDy);
                        % else
                        %     error('Number of voxels in y is not an integer.')
                        % end
                        % if abs(round(nVPDz) - nVPDz) < 1e-12
                        %     nVPDz = round(nVPDz);
                        % else
                        %     error('Number of voxels in z is not an integer.')
                        % end
                        % nVoxels = double(uint32(nVPDx) * uint32(nVPDy) * uint32(nVPDz));
                        %
                        % %shape_cropped = [nVPDx, nVPDy, nVPDz];
                        %
                        % bead_struct.Beads = cropToDomain(bead_struct.Beads, ...
                        %                                     nVoxelsOld, crop_mask);
                        % %%%%% THERE'S A BUG HERE... NEED TO RE-ORGANIZE
                        % bead_struct.EdgeIndices = cropToDomain(bead_struct.EdgeIndices, ...
                        %                                         nVoxelsOld, crop_mask);
                        % bead_struct.Shell = cropToDomain(bead_struct.Shell, ...
                        %                                     nVoxelsOld, crop_mask);
                        %
                        % bead_struct.Beads = bead_struct.Beads(~cellfun(@isempty, bead_struct.Beads));
                        % bead_struct.EdgeIndices = bead_struct.EdgeIndices(~cellfun(@isempty, bead_struct.EdgeIndices));
                        % bead_struct.Shell = bead_struct.Shell(~cellfun(@isempty, bead_struct.Shell));
                        %
                        % shape = [nVPDx, nVPDy, nVPDz];
                        % % Update bead struct
                        % bead_struct.Shape     = shape;
                        % bead_struct.DomainMin = [domain(1), domain(3), domain(5)];
                        % bead_struct.AllBeads  = cell2mat(bead_struct.Beads);
                        % bead_struct.AllEdges  = cell2mat(bead_struct.EdgeIndices);
                        % bead_struct.AllShells = cell2mat(bead_struct.Shell);

                    end
                    % Compute particle diameter distribution
                    diameters = particleDiam(bead_struct.Beads, voxels, dx);

                    % Storing data
                    data.beads        = bead_struct.Beads;
                    data.edgeIndices  = bead_struct.EdgeIndices;
                    data.shell        = bead_struct.Shell;
                    data.allBeads     = bead_struct.AllBeads;
                    data.allEdges     = bead_struct.AllEdges;

                    bead_centers = bead_data(:, 1:3);
                    num_beads = length(bead_struct.Beads);
                    bead_diams = diameters;

                    % fprintf('%30s %1.1e\n', 'Number of voxels:', nVoxels);
                    % fprintf(runtimes_file, '%30s %1.1e\n', 'Number of voxels:', nVoxels);
                    % fprintf('%30s %i\n', 'Number of beads:', num_beads);
                    % fprintf(runtimes_file, '%30s %i\n', 'Number of beads:', num_beads);
                    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Label bead domain:');

                case 'labeled'
                    beads = dlmread(domain_file, ',', row_cntr-1, 0);

                    % Create 3D grid of voxel coordinates (center of voxel cube)
                    voxels = centeredGrid3D(domain, dx);

                    % Shape
                    nVPDx = (domain(2) - domain(1)) / dx;
                    nVPDy = (domain(4) - domain(3)) / dx;
                    nVPDz = (domain(6) - domain(5)) / dx;

                    if abs(round(nVPDx) - nVPDx) < 1e-12
                        nVPDx = round(nVPDx);
                    else
                        error('Number of voxels in x is not an integer.')
                    end
                    if abs(round(nVPDy) - nVPDy) < 1e-12
                        nVPDy = round(nVPDy);
                    else
                        error('Number of voxels in y is not an integer.')
                    end
                    if abs(round(nVPDz) - nVPDz) < 1e-12
                        nVPDz = round(nVPDz);
                    else
                        error('Number of voxels in z is not an integer.')
                    end

                    shape = [nVPDx, nVPDy, nVPDz];
                    nVoxels = nVPDx * nVPDy * nVPDz;

                    if crop_percent < 1
                        [voxels, domain, voxel_indices] = cropBeadDomain_voxels(voxels, domain, ...
                                                                 crop_percent, dx);
                        % Shape
                        nVPDx = (domain(2) - domain(1)) / dx;
                        nVPDy = (domain(4) - domain(3)) / dx;
                        nVPDz = (domain(6) - domain(5)) / dx;

                        if abs(round(nVPDx) - nVPDx) < 1e-12
                            nVPDx = round(nVPDx);
                        else
                            error('Number of voxels in x is not an integer.')
                        end
                        if abs(round(nVPDy) - nVPDy) < 1e-12
                            nVPDy = round(nVPDy);
                        else
                            error('Number of voxels in y is not an integer.')
                        end
                        if abs(round(nVPDz) - nVPDz) < 1e-12
                            nVPDz = round(nVPDz);
                        else
                            error('Number of voxels in z is not an integer.')
                        end
                        shape_cropped = [nVPDx, nVPDy, nVPDz];
                        nVoxels = nVPDx * nVPDy * nVPDz;

                        % collect rows in cropped domain
                        [beads_sorted, sort_ind] = sort(beads(:, 1));
                        matching_rows = fastIntersect(beads_sorted, voxel_indices, ...
                                                      'bool vec');
                        beads = sortrows(beads(sort_ind(matching_rows), :), 2);
                        % convert row indices in 1st column of beads to cropped indices
                        [xf, yf, zf] = ind2sub(shape, (beads(:, 1) - min(beads(:, 1)) + 1));
                        beads = [all_sub2ind(shape_cropped, [xf yf zf]), beads(:, 2)];

                        shape = shape_cropped;
                    end

                    data.numVoxels = nVoxels;
                    data.dx = dx;
                    data.domain = domain;
                    data.shape = shape;
                    % fprintf('%30s %1.1e\n', 'Number of voxels:', nVoxels);
                    % fprintf(runtimes_file, '%30s %1.1e\n', 'Number of voxels:', nVoxels);

                    bead_struct = struct;
                    % Save number of beads
                    num_beads = length(find(diff(beads(:,2)) ~= 0)) + 1;
                    % fprintf('%30s %i\n', 'Number of beads:', num_beads);
                    % fprintf(runtimes_file, '%30s %i\n', 'Number of beads:', num_beads);

                    % Additional information needed to store full bead data
                    domMin     = voxels(1, :);
                    shell_voxs = shell_thickness / dx;
                    bead_struct.Beads       = cell(num_beads, 1);
                    bead_struct.EdgeIndices = cell(num_beads, 1);
                    bead_struct.Shell       = cell(num_beads, 1);

        %             % Sort (data should already be sorted)
        %             beads = sortrows(beads, 2);

                    % gather indices of each bead
                    iter = 1;
                    bead_current = beads(1, 2);
                    bead_num     = 1;
                    bead_centers = zeros(num_beads, 3);
                    while bead_num <= num_beads
                        % collect indices of bead number 'bead_current'
                        beg_ind = iter;
                        while iter < size(beads, 1) && beads(iter + 1, 2) == bead_current
                            iter = iter + 1;
                        end
                        end_ind  = iter;
                        row_inds = beg_ind : end_ind;
                        bead_struct.Beads{bead_num} = [bead_struct.Beads{bead_num}; beads(row_inds, 1)];
                        % Find bead centers to compute convex hull
                        bead_centers(bead_num, :) = mean(voxels(bead_struct.Beads{bead_num}, :), 1);

                        % Additional bead information
                        % Store edge voxels using bwdist
                        % shift domain voxels to (1, 1, 1);
                        % divide by dx to get matrix voxel as opposed to real-space point
                        bead_voxs = voxels(bead_struct.Beads{bead_num}, :);
                        bead_vox0 = floor((bead_voxs - domMin) / dx) + 1;
                        % add 1 voxel buffer
                        xmin_i = min(bead_vox0(:, 1)) - 1;
                        ymin_i = min(bead_vox0(:, 2)) - 1;
                        zmin_i = min(bead_vox0(:, 3)) - 1;
                        xmax_i = max(bead_vox0(:, 1)) + 1;
                        ymax_i = max(bead_vox0(:, 2)) + 1;
                        zmax_i = max(bead_vox0(:, 3)) + 1;

                        xIdx   = xmin_i : xmax_i;
                        yIdx   = ymin_i : ymax_i;
                        zIdx   = zmin_i : zmax_i;
                        shape_box = [length(xIdx), length(yIdx), length(zIdx)];

                        [xxIdx, yyIdx, zzIdx] = ndgrid(xIdx, yIdx, zIdx);
                        allIdx = [xxIdx(:) yyIdx(:) zzIdx(:)];

                        allLinIdx         = all_sub2ind(shape, allIdx);
                        bead_vox0_bb      = bead_vox0 - allIdx(1, :) + 1; % bounding box
                        bead_boxind       = all_sub2ind(shape_box, bead_vox0_bb);

                        % Store edge voxels using bwdist
                        bin_imbead = ones(shape_box);
                        bin_imbead(bead_boxind) = 0;
                        edt_edgebead = dx * bwdist(bin_imbead);
                        edt_edge_indbead = find(edt_edgebead <= (dx * sqrt(2)) + 1e-5);
                        edge_indbead = intersect(edt_edge_indbead, bead_boxind);
                        bead_struct.EdgeIndices{bead_num} = allLinIdx(edge_indbead);

                        % Store a shell cloud of points that is shell_thickness thick
                        edt_shell_ind = find(edt_edgebead <= (shell_voxs * dx * sqrt(2)) + 1e-5);
                        shell_ind = intersect(edt_shell_ind, bead_boxind);
                        bead_struct.Shell{bead_num} = allLinIdx(shell_ind);

                        % Update counters
                        iter = iter + 1;
                        bead_num = bead_num + 1;
                        if iter <= size(beads, 1)
                            bead_current = beads(iter, 2);
                        else
                            bead_num = num_beads + 1;
                        end
                    end
                    % Additional information for struct
                    bead_struct.Shape     = shape;
                    bead_struct.DomainMin = domMin;

                    % Compute particle diameter distribution
                    diameters = particleDiam(bead_struct.Beads, voxels, dx);

                    % Get cumulative bead indices list
                    bead_struct.AllBeads = [];
                    bead_struct.AllEdges = [];
                    for j = 1 : num_beads
                        bead_struct.AllBeads = [bead_struct.AllBeads; bead_struct.Beads{j}];
                        bead_struct.AllEdges = [bead_struct.AllEdges; bead_struct.EdgeIndices{j}];
                    end
                    bead_struct.AllEdges = unique(bead_struct.AllEdges);
                    % Store
                    data.beads        = bead_struct.Beads;
                    data.edgeIndices  = bead_struct.EdgeIndices;
                    data.shell        = bead_struct.Shell;
                    data.allBeads     = bead_struct.AllBeads;
                    data.allEdges     = bead_struct.AllEdges;

                    bead_diams = diameters;
                otherwise
                    error('Data file does not have a known data type.')
            end
        otherwise
            error('Unexpected file extension: %s', fExt);
    end

    tElapsed = toc(tStart);
    time_log(timeLogIdx).Name = 'Process input data';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % This is used later to look up what bead owns what voxels.
    voxelToBeadMap = VoxelToBeadMap(data.beads);

    %%%%%%%*************************************************************%%%%%%%
    %%%%%********************** CROP TO CONVEX HULL **********************%%%%%
    %%%%%%%*************************************************************%%%%%%%

    tStart = tic;

    K = convhull(bead_centers);
    in_log = inhull(voxels, bead_centers, K);
    in_log3D = reshape(in_log, shape);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Crop to convex hull:');

    time_log(timeLogIdx).Name = 'Crop to convex hull';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    bead_voxels = false(nVoxels, 1);

    % Obtain void space voxels and exterior voxels (i.e., non-bead,
    % non-void space voxels)
    vsVoxel_log = in_log;
    vsVoxel_log(data.allBeads) = false;
    bead_voxels(data.allBeads) = true;

    exterior_voxels = (1 : nVoxels)';
    exterior_voxels = exterior_voxels(~vsVoxel_log & ~bead_voxels);
    domainBoundary_voxels = wallIndices(voxels, dx);

    vsVoxels = (1 : nVoxels)';
    vsVoxels = vsVoxels(vsVoxel_log);
    num_VsVoxels = length(vsVoxels);

    data.voidVoxels = vsVoxels;

    % Determine edge voxels of void space
    edge_ind = edgeVoidVoxels(bead_voxels, vsVoxels, shape, dx, 1);

    % Calculate void volume fraction within void space
    porosity = num_VsVoxels / sum(in_log);

    dx
    nVoxels

    % Calculate void area fraction along 2D slice
    zgMin = ceil(data.domain(5) / dx);
    zgMax = ceil(data.domain(6) / dx);
    avoid_edges = ceil((zgMax - zgMin) / 6);
    slices = floor(linspace(zgMin + avoid_edges, zgMax - avoid_edges, num_2D_slices));
    if slices(2) - slices(1) < dx
        slices = zgMin + avoid_edges : dx : zgMax - avoid_edges;
    end

%         stats2Dslice = cell(numel(slices), 1);
%         num_voxs2Dslice = shape(1) * shape(2);
%         shape2Dslice = [shape(1), shape(2)];
    slice_vals = zeros(numel(slices), 1);
    for i = 1 : numel(slices)
        % Void area
        z_plane_log     = ceil(voxels(:, 3) ./ dx) == slices(i);
        vs_plane_nVoxs  = vsVoxel_log & z_plane_log; % void space only
        tot_plane_nVoxs = in_log & z_plane_log; % complete convex hull
        slice_vals(i)   = sum(vs_plane_nVoxs) / sum(tot_plane_nVoxs);
%
%             % 2D subunits
%             bin_image2D = ones(shape2Dslice);
%             [a, b, ~]   = ind2sub(shape, find(vs_plane_nVoxs));
%             vs_2Dslice = sub2ind([shape(1), shape(2)], a, b);
%             vs_2Dslicelog = false(shape2Dslice);
%             vs_2Dslicelog(vs_2Dslice) = true;
%             bin_image2D(vs_2Dslicelog) = 0;
%             edt_2Dslice = dx * bwdist(bin_image2D);
%             rooms2Dslice = (edt_2Dslice > (mean(bead_diams) * (hall_cutoff / 100))); % hall_cutoff
%
%             if sum(rooms2Dslice(:)) > 0
%                 % isolate islands of 'rooms'
%                 conncomp_2Dslice = bwconncomp(rooms2Dslice);
%                 room_islands2Dslice = conncomp_2Dslice.PixelIdxList;
%                 rm_sizes = cell2mat(cellfun(@numel, room_islands2Dslice, 'UniformOutput', false));
%                 rmv_rms  = rm_sizes <= 1;
%                 room_islands2Dslice(rmv_rms) = [];
%                 num_islands2Dslice  = length(room_islands2Dslice);
%                 room2Dslice_key     = zeros(num_voxs2Dslice, 1);
%                 % Associate hallway voxels to corresponding rooms to form pore subunits
%                 % create binary image of islands
%                 nn2D_binImage = zeros(shape2Dslice);
%                 for j = 1 : num_islands2Dslice
%                     nn2D_binImage(room_islands2Dslice{j})   = 1;
%                     room2Dslice_key(room_islands2Dslice{j}) = j;
%                 end
%                 % use nearest-index info to associate voxel to nearest island
%                 [~, nn2D_idx] = bwdist(nn2D_binImage);
%                 % remove non-void space voxels
%                 nn2D_idx = nn2D_idx(:);
%                 nn2D_idx(~vs_2Dslicelog) = [];
%
%                 % form subunits
%                 nn2D_idx_subs = room2Dslice_key(nn2D_idx);
%
%                 subs_byHall2Dslice = room_islands2Dslice';
%                 for j = 1 : num_islands2Dslice
%                     subs_byHall2Dslice{j} = sort([subs_byHall2Dslice{j}; vs_2Dslice(nn2D_idx_subs == j)]);
%                 end
% %                 % for each sub, grab peaks
% %                 gaus  = imgaussfilt(edt_2Dslice);
% %                 gaus2 = imgaussfilt(gaus);
% %                 pks   = imregionalmax(gaus2);
%
%                 % reshape data into labeled domain for regionprops fnc
%                 labeled_input = zeros(shape2Dslice);
%                 for j = 1 : num_islands2Dslice
%                     labeled_input(subs_byHall2Dslice{j}) = j;
%                 end
%                 anseth_meth = regionprops(labeled_input, 'Area', 'MajorAxisLength', 'MinorAxisLength');
%                 num_am      = numel(anseth_meth);
%                 am_area     = arrayfun(@(x) anseth_meth(x).Area, 1 : num_am)';
%                 am_diam     = sqrt((am_area / pi)) * 2;
%                 am_maxaxis  = arrayfun(@(x) anseth_meth(x).MajorAxisLength, 1 : num_am)';
%                 am_minaxis  = arrayfun(@(x) anseth_meth(x).MinorAxisLength, 1 : num_am)';
%                 stats2Dslice{i} = [am_area, am_diam, am_maxaxis, am_minaxis];
%             else
%                 stats2Dslice{i} = [];
%             end
    end
%         stats2Dslice_all = cell2mat(stats2Dslice);
%         % write to excel
%         if strcmp(data_filename(1:8), 'beadInfo') || strcmp(data_filename(1:13), 'labeledDomain')
%         excel2Dslice_filename = replace(data_filename, {'beadInfo_', 'labeledDomain_', '.dat', '.txt', '.json'}, ...
%                                                     {'stats2Dslice_', 'stats2Dslice_', '.xlsx', '.xlsx', '.xlsx'});
%         else
%             excel2Dslice_filename = replace(data_filename, {'.dat', '.txt', '.json'}, ...
%                                                     {'.xlsx', '.xlsx', '.xlsx'});
%             excel2Dslice_filename = horzcat('stats2Dslice_', excel2Dslice_filename);
%         end
%         excel2Dslice_filename = ['./outputs/', excel2Dslice_filename];
%         writecell({'Area', 'Circle Diameter', 'Major Axis', 'Minor Axis'}, excel2Dslice_filename);
%         writematrix(stats2Dslice_all, excel2Dslice_filename, 'WriteMode', 'append');
%
%         return;

    porosity2D = mean(slice_vals);

    void_vol_fract     = porosity;
    particle_fract     = 1 - porosity;
    void_area_fract    = porosity2D;

    %%%%%%%*************************************************************%%%%%%%
    %%%%%***************** SUBUNITS BY HALLWAY THRESHOLD *****************%%%%%
    %%%%%%%*************************************************************%%%%%%%

    % Form subunits using hall_cutoff threshold
    tStart = tic;
    bin_image = zeros(shape);
    bin_image(data.allBeads) = 1;

    % compute EDT on binary image
    edt = dx * bwdist(bin_image);
    % filter by EDT threshold that removes 'narrow hallways'
    rooms = (edt > hall_cutoff) & in_log3D;

    hallways = (edt <= hall_cutoff) & reshape(vsVoxel_log, shape);
    data.hallways_indices = find(hallways);

    if sum(rooms(:)) > 0
        % isolate islands of 'rooms'
        room_islands = bwconncomp(rooms);
        num_islands  = length(room_islands.PixelIdxList);
        room_key     = zeros(nVoxels, 1);
        % Associate hallway voxels to corresponding rooms to form pore subunits
        % create binary image of islands
        nn_binImage = zeros(shape);
        for i = 1 : num_islands
            nn_binImage(room_islands.PixelIdxList{i}) = 1;
            room_key(room_islands.PixelIdxList{i}) = i;
        end
        % use nearest-index info to associate voxel to nearest island
        [~, nn_idx] = bwdist(nn_binImage);
        % remove non-void space voxels
        nn_idx = nn_idx(:);
        nn_idx(~vsVoxel_log) = [];

        % form subunits
        nn_idx_subs = room_key(nn_idx);

        subs_byHall = room_islands.PixelIdxList';
        for i = 1 : num_islands
            subs_byHall{i} = sort([subs_byHall{i}; vsVoxels(nn_idx_subs == i)]);
        end
    else
        subs_byHall = {};
    end

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Form subunits via hall cutoff:');

    time_log(timeLogIdx).Name = 'Form subunits via hall cutoff';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    %%%%%%%*******************************************************%%%%%%%
    %%%%%***************** NEAREST BEADS PER VOXEL *****************%%%%%
    %%%%%%%*******************************************************%%%%%%%

    tStart = tic;

    [nearestBeads, sparseBeadEDTs, data.EDT] = ...
        findNearestBeads(shape, data.dx, data.beads, vsVoxels, sqrt(2));
    data.EDT(~vsVoxel_log) = -1;

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Find nearest beads:');

    time_log(timeLogIdx).Name = 'Find nearest beads';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % Voxel Type
    tStart = tic;
    % Obtaining peaks, ridges1D, and ridges2D data
    % e.g., .indices  - column array containing indices of voxels
    %       .beads    - matrix where columns refer to unique peaks/ridges
    %                   and rows are nearest beads in order
    %       .numBeads - column array containing the # of beads that are
    %                   equally as close
    [medax, ridges2D] = determineVoxelType(nearestBeads, vsVoxels);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Determine voxel type:');

    time_log(timeLogIdx).Name = 'Determine voxel type';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    tStart = tic;
    % Determine regions of available void space depending on cell/species size
    % The 4 cell/species diameter categories are: 0, 10, 30, 60 um
    % The values in these categories are a set constant for consistency/comparison
    % among data. They were chosen by calibrating for a 40 um and 100 um
    % rigid homogeneous domain. Region bin-sizes to store data were also
    % calibrated for the same domains.
    % Bin volume cutoffs are in volume (pL = um^3 / 1000)
    bin_volume_cutoffs = [0, 100, 1000]; % in picoliters
    species_cutoffs    = [0, 10, 30, 60]; % in um

    % for each species size cutoff
    % region_bins key:
    %   small bin  = 1
    %   medium bin = 2
    %   large bin  = 3
    %   empty bin  = 4
    region_vols     = cell(numel(species_cutoffs), 1);
    region_vols_ind = cell(numel(species_cutoffs), 1);
    region_bins     = cell(numel(species_cutoffs), 1);
    region_data     = cell(numel(species_cutoffs), 1);
    regions_key_output = cell(numel(species_cutoffs), 1);
    % do not include these indices in analysis because not connected to
    % main void space backbone
    throw_away_bb = [];
    num_regions_removed = 0;
    for i = 1 : numel(species_cutoffs)
        % for each species size, determine the medial axes voxels with
        % EDTs above size / 2 (i.e., radius)
        medax_log = false(shape);
        medax_log(medax.indices) = true;
        these_regions     = (edt >= (species_cutoffs(i) / 2)) & medax_log;
        remaining_regions = (edt < (species_cutoffs(i) / 2)) & medax_log;
        clear('medax_log')
        % locate unique islands
        region_islands = bwconncomp(these_regions);
        tot_islands    = numel(region_islands.PixelIdxList);
        if isempty(region_islands.PixelIdxList)
            region_vols{i} = 0;
            region_bins{i} = 4;
            regions_key_output{i} = {};
        elseif species_cutoffs(i) == 0
            s0 = i;
            % honestly i think this is over-kill and unnecessary...but i guess we're doing this
            % We're removing regions of the void space that are not connected to the main backbone of the void space
            % normalize data by species_cutoff = 0 so that any regions
            % not part of the largest region are not considered
            [~, which_ind] = max(cellfun(@(x) numel(x), region_islands.PixelIdxList));
            % store indices not part of the largest region
            which_inds_all = [1 : which_ind - 1, which_ind + 1 : tot_islands];
            keep_log = true(nVoxels, 1);
            if ~isempty(which_inds_all)
                throw_away_bb = sort(vertcat(region_islands.PixelIdxList{which_inds_all}));
                for j = 1 : numel(which_inds_all)
                    keep_log(regionVolume(region_islands.PixelIdxList{which_inds_all(j)}, ...
                                    data.EDT, shape, dx)) = false;
                end
                num_regions_removed = tot_islands - numel(which_inds_all);
            else
                num_regions_removed = 0;
            end
            % save volume indices
            region_vols_ind{i} = (1 : nVoxels)'; % all void space voxels
            region_vols_ind{i} = region_vols_ind{i}(keep_log & vsVoxel_log);
            % save volume
            region_vols{i} = round2sigdig(numel(region_vols_ind{i}) * dx^3, dx) / 1000;
            region_bins{i} = 3;
            % save regions key
            regions_key_output{i}{1} = find(keep_log);
        else
            if tot_islands == num_regions_removed + 1
                % No change from species_cutoff = 0
                % save volume
                region_vols{i} = region_vols{s0};
                region_bins{i} = 3;
                regions_key_output{i} = zeros(nVoxels, 1);
                regions_key_output{i}(keep_log) = 1;
            else
                % remove regions of the void space that are not connected to the main backbone of the void space
                remove_regions = false(tot_islands, 1);
                for j = 1 : tot_islands
                    if fastIntersect(region_islands.PixelIdxList{j}, throw_away_bb, 'bool')
                        remove_regions(j) = true;
                    end
                end
                region_islands.PixelIdxList(remove_regions) = [];
                tot_islands = tot_islands - sum(remove_regions);

                % Now find volume of each region
                region_key = zeros(nVoxels, 1);
                % Associate void voxels to corresponding regions to create local volume of islands
                % Create binary image of islands
                region_binImage = zeros(shape);
                region_binImage(remaining_regions) = 1;
                region_key(remaining_regions)      = -1;
                for j = 1 : tot_islands
                    region_binImage(region_islands.PixelIdxList{j}) = 1;
                    region_key(region_islands.PixelIdxList{j})      = j;
                end

                % extracting volume
                region_vols_ind{i} = cell(tot_islands, 1);
                region_vols{i}     = zeros(tot_islands, 1);
                region_bins{i}     = zeros(tot_islands, 1);
                for j = 1 : tot_islands
                    region_vols_ind{i}{j} = regionVolume(region_islands.PixelIdxList{j}, data.EDT, shape, dx);
                    r_vol = round2sigdig(numel(region_vols_ind{i}{j}) * dx^3, dx) / 1000;
                    if r_vol < bin_volume_cutoffs(2)
                        region_bins{i}(j)  = 1;
                    elseif r_vol < bin_volume_cutoffs(3)
                        region_bins{i}(j)  = 2;
                    else
                        region_bins{i}(j) = 3;
                    end
                    region_vols{i}(j) = r_vol;
                    regions_key_output{i}{j} = region_vols_ind{i}{j};
                end
%                     region_vols{i}(end) = region_vols{s0} - sum(region_vols{i});
%                     region_bins{i}(end) = 4;
            end
        end
        num_regions = size(region_vols{i}, 1);
        region_data{i} = [(1 : num_regions)', region_vols{i}(:), region_bins{i}(:)];
    end
    clear('region_binImage')

    data.Regions.volume_data = region_data;
    data.Regions.regions_key = regions_key_output;

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Regions based on species size:');

    time_log(timeLogIdx).Name = 'Regions based on species size';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % get edge voxels of beads at edge of void space
    bead_struct.EdgeBeadIndices = edgeBeadVoidVoxels(in_log, data.edgeIndices, shape, dx, 1);

    %%%%%%%*****************************************************%%%%%%%
    %%%%%***************** FINDING PEAKS ON EDGE *****************%%%%%
    %%%%%%%*****************************************************%%%%%%%

    % Edge (of Void Space) Voxels
    tStart = tic;

    [nearestBeadsEdge, sparseBeadEdgeEDTs, data.edgeEDT] = ...
        findNearestBeads(shape, data.dx, bead_struct.EdgeBeadIndices, edge_ind, sqrt(2));
    data.edgeEDT = data.edgeEDT(edge_ind);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Find nearest beads edge:');

    time_log(timeLogIdx).Name = 'Find nearest beads edge';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % Voxel Type
    tStart = tic;
    % Create ridges1D
    peaks_e    = struct;
    % Obtaining peaks, ridges1D, and ridges2D data
    % e.g., .indices  - column array containing indices of voxels
    %       .beads    - matrix where columns refer to unique peaks/ridges
    %                   and rows are nearest beads in order
    %       .numBeads - column array containing the # of beads that are
    %                   equally as close
    [medax_e, ridges_e] = determineVoxelType(nearestBeadsEdge, edge_ind);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Determine voxel type edge:');

    time_log(timeLogIdx).Name = 'Determine voxel type edge';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % Construct the full EDT
    edt_full_e_key = zeros(nVoxels, 1);
    edt_full_e_key(edge_ind) = data.edgeEDT;

    % Find peaks on edges
    edge_bw = zeros(shape);
    edge_bw(medax_e.indices) = 1;
    peak_clusters_e = bwconncomp(edge_bw);
    peaks_e.num = length(peak_clusters_e.PixelIdxList);
    % Find highest EDT in each cluster
    peaks_e.indices = zeros(peaks_e.num, 1);
    for i = 1 : peaks_e.num
        [~, ind] = max(edt_full_e_key(peak_clusters_e.PixelIdxList{i}));
        peaks_e.indices(i) = peak_clusters_e.PixelIdxList{i}(ind);
    end

    % Apply hallway cutoff to form island clusters
    roomsEdge_2D = ridges_e.cum_indices(edt_full_e_key(ridges_e.cum_indices) > hall_cutoff);
    if sum(roomsEdge_2D(:)) > 0
        edge_bw = zeros(shape);
        edge_bw(roomsEdge_2D) = 1;
        edge_bw(medax_e.indices) = 1;
    else
        edge_bw = zeros(shape);
        edge_bw(peaks_e.indices) = 1;
    end

    % isolate islands of 2D 'rooms' (or subs)
    room_islands_2D = bwconncomp(edge_bw);
    edgesubs2comb_2Dind = room_islands_2D.PixelIdxList;

    %rm_sizes    = cell2mat(cellfun(@numel, edgesubs2comb_2Dind, 'UniformOutput', false));
    num_2Dedgesubs = numel(edgesubs2comb_2Dind);
    edge_2Dsubs    = cell(num_2Dedgesubs, 1);

    subs2D_key   = zeros(num_2Dedgesubs, 1);
    % Associate hallway voxels to corresponding rooms to form edge subunit backbones
    % Create binary image of islands
    nnedge2D_binIm = zeros(shape);
    for j = 1 : num_2Dedgesubs
        nnedge2D_binIm(edgesubs2comb_2Dind{j}) = 1;
        subs2D_key(edgesubs2comb_2Dind{j})     = j;
    end
    % use nearest-index info to associate voxel to nearest island
    [~, nn2D_ind] = bwdist(nnedge2D_binIm); % jankily using distance through 3D space even though we're living on the 2D surface...
    % remove non-edge voxels
    edge_log = false(nVoxels, 1);
    edge_log(edge_ind) = true;
    nn2D_ind = nn2D_ind(:);
    nn2D_ind(~edge_log) = [];

    % form subunits
    nn2D_ind_subs = subs2D_key(nn2D_ind);

    for j = 1 : num_2Dedgesubs
        edge_2Dsubs{j} = struct;
        edge_2Dsubs{j}.indices = edgesubs2comb_2Dind{j};
        edge_2Dsubs{j}.indices = sort([edge_2Dsubs{j}.indices; edge_ind(nn2D_ind_subs == j)]);
        edge_2Dsubs{j}.ridges  = edgesubs2comb_2Dind{j};
    end

    % Compute accessible ligand (in micromoles) of edge subs by identifying edge
    % voxels of neighbor beads and multiplying by thickness of shell
    beads_by_num = zeros(shape);
    for i = 1 : numel(bead_struct.EdgeBeadIndices)
        beads_by_num(bead_struct.EdgeBeadIndices{i}) = i;
    end

    beads_by_num_padded = padarray(beads_by_num, [1 1 1], 0, 'both');
    for i = 1 : num_2Dedgesubs
        num_2Dedgevoxs = numel(edge_2Dsubs{i}.indices);
        edge_2Dshellbool = false(num_2Dedgevoxs, 1);
        edge_2Dshellinds = zeros(num_2Dedgevoxs, 1);
        for j = 1 : num_2Dedgevoxs
            bounce = false;
            % check if edge voxel neighbors a bead
            [x,y,z] = ind2sub(shape, edge_2Dsubs{i}.indices(j));
            % shift indices to adjust for pad
            x = x + 1;
            y = y + 1;
            z = z + 1;
            for k = x - 1 : x + 1
                for l = y - 1 : y + 1
                    for m = z - 1 : z + 1
                        if beads_by_num_padded(k,l,m) ~= 0
                            edge_2Dshellbool(j) = true;
                            % edge_beadnum(j)   = beads_by_number(k,l,m);
                            edge_2Dshellinds(j) = sub2ind(shape, k - 1, l - 1, m - 1);
                            bounce = true;
                        end
                        if bounce
                            break;
                        end
                    end
                    if bounce
                        break;
                    end
                end
                if bounce
                    break;
                end
            end
        end
        edge_2Dshellinds = edge_2Dshellinds(edge_2Dshellinds ~= 0);
        edge_2Dsubs{i}.touchingBeadindices = edge_2Dshellinds;
        % Store edge ligand information
        num_2Dedgeshellvoxs = sum(edge_2Dshellbool);
        SA_2Dedgeshell = num_2Dedgeshellvoxs * dx^2;
        tot_equiv_vol2D = SA_2Dedgeshell * shell_thickness; % in um^3
        % using conversion:
        % 1 um^3 = 1e-15 L
        RGD_conc = 500;
        edge_2Dsubs{i}.RGDaccessible = RGD_conc * tot_equiv_vol2D * 1e-15;
    end

    % Match up peaks to edge subunit clusters
    edgesubs2comb_cum = cell2mat(arrayfun(@(x) [edgesubs2comb_2Dind{x}, ...
                                          x * ones(length(edgesubs2comb_2Dind{x}), 1)]', ...
                                          1 : length(edgesubs2comb_2Dind), ...
                                          'UniformOutput', false))';

    knnsearch_ptCloud_e = knnsearch(voxels(edgesubs2comb_cum(:, 1), :), ...
                                    voxels(peaks_e.indices, :), ...
                                    'k', 1, 'NSMethod', 'kdtree');

    % convert knnsearch output to corresponding edge island number
    peaks_e.subs_key = edgesubs2comb_cum(knnsearch_ptCloud_e, 2);

%         else
%             edgesubs2comb_2Dind = arrayfun(@(x) {x}, peaks_e.indices);
%             peaks_e.subs_key    = (1 : length(peaks_e.indices))';
%         end
%         else
%             edge_2Dsubs         = {};
%             edgesubs2comb_2Dind = {};
%             num_2Dedgesubs      = 0;
%         end

    data.SurfaceSubs = edge_2Dsubs;

    % Assign peak to edge_2Dsubs
    for i = 1 : num_2Dedgesubs
        edge_2Dsubs{i}; % NEED TO FINISH
    end

    % Create 'door' for each peak
    % per peak, find voxel on nearest 3 beads; each column is a
    % different bead
    peaks_e.bead_voxs = zeros(peaks_e.num, 3);
    peaks_e.center    = cell(peaks_e.num, 1);
    peaks_e.radius    = cell(peaks_e.num, 1);
    peaks_e.v1        = cell(peaks_e.num, 1);
    peaks_e.v2        = cell(peaks_e.num, 1);
    peaks_e.nvec      = cell(peaks_e.num, 1);
    for i = 1 : peaks_e.num
        ind = binarySearch(medax_e.indices, peaks_e.indices(i));
        allNearBeads = getNearestBeads(medax_e.beads, ind);
        allDists = zeros(size(allNearBeads));
        for ii = 1 : numel(allDists)
            allDists(ii) = sparseBeadEdgeEDTs{allNearBeads(ii)}(peaks_e.indices(i));
        end
        [~, sortPerm] = sortrows([allDists, allNearBeads]);
        beads_current = allNearBeads(sortPerm(1 : 3));

        % checking nearest 3 beads since need 3 pts to form a plane
        for h = 1 : 3
            edgebead_voxs = bead_struct.EdgeBeadIndices{beads_current(h)};
            norm_vec = zeros(length(edgebead_voxs), 1);
            for j = 1 : length(edgebead_voxs)
                norm_vec(j) = norm(voxels(edgebead_voxs(j), :) - voxels(peaks_e.indices(i), :));
            end
            [~, minDist_idx] = min(norm_vec);
            peaks_e.bead_voxs(i, h) = edgebead_voxs(minDist_idx);
        end
        [center, radius, v1, v2] = circlefit3d(voxels(peaks_e.bead_voxs(i, 1), :), ...
                                               voxels(peaks_e.bead_voxs(i, 2), :), ...
                                               voxels(peaks_e.bead_voxs(i, 3), :));
        peaks_e.center{i} = center;
        peaks_e.radius{i} = radius;
        peaks_e.v1{i}     = v1;
        peaks_e.v2{i}     = v2;
        if ~isempty(v1 + v2)
            peaks_e.nvec{i} = cross(v1, v2);
        else
            peaks_e.nvec{i} = [NaN; NaN; NaN];
        end
    end

%         % For each cluster, convert indices (based on shape) to indices
%         % based on edge_ind. Grab nearest beads for each voxel in
%         % cluster and compare to make sure they share the same nearest beads
%         beads_e = cellfun(@(x) medax_e.beads(:, fastIntersect(medax_e.indices, x, 'indices')), ...
%                           peak_clusters_e.PixelIdxList, 'UniformOutput', false);
%         peaks_e.indices = zeros(num_peaks_e, 1);
%         for i = 1 : length(beads_e)
%             if size(unique(beads_e{i}), 2) == size(beads_e{i}, 1)
%                 % all peaks in the cluster share the same bead nei
%                 [~, ind] = max(edt_full_e_key(peak_clusters_e.PixelIdxList{i}));
%                 peaks_e.indices(i) = peak_clusters_e.PixelIdxList{i}(ind);
%             else
%
%             end
%         end

    %%%%%%%**********************************************************%%%%%%%
    %%%%%***************** IDENTIFYING RIDGES & PEAKS *****************%%%%%
    %%%%%%%**********************************************************%%%%%%%

    tStart = tic;
    % Find unique 2D ridges based on common neighboring beads
    % e.g., .indices  - cell array where each cell contains the indices of
    %                   voxels that make up that unique ridge
    %       .beads    - cell array containing the 2 nearest beads in order

    [c, ~, ic] = unique(ridges2D.cum_beads', 'rows');
    ridges2D.indices = accumarray(ic, ridges2D.cum_indices, [], @(r) {sort(r)});
    ridges2D.beads = mat2cell(c, ones(size(c, 1), 1));
    ridges2D.beadNeighbors = vertcat(ridges2D.beads{:});

    % Obtain crawl space info
    num_r2D = numel(ridges2D.indices);
    crawl_width = zeros(num_r2D, 1);
    crawl_inds = cell(num_r2D, 1);
    touching_beads_log = false(num_r2D, 1);
    for i = 1 : num_r2D
        val = min(edt(ridges2D.indices{i}));
        % assume [val = dx] means touching particles
        if val <= dx
            crawl_width(i) = 0;
            touching_beads_log(i) = true;
        else
            crawl_width(i) = val * 2;
        end
        crawl_inds{i} = ridges2D.indices{i}(edt(ridges2D.indices{i}) == val);
    end
    ridges2D.crawl_inds = crawl_inds;
    ridges2D.crawl_widths = crawl_width;

    % Create sparse full adjacency matrix
    beads_adjacency = sparse(ridges2D.beadNeighbors(:, 1), ...
                             ridges2D.beadNeighbors(:, 2), ...
                             ones(size(ridges2D.beadNeighbors, 1), 1), ...
                             num_beads, num_beads);
    % Create sparse adjacency matrix of only touching particles
    only_touching = ridges2D.beadNeighbors(touching_beads_log, :);
    touching_adjacency = sparse(only_touching(:, 1), ...
                                only_touching(:, 2), ...
                                ones(size(only_touching, 1), 1), ...
                                num_beads, num_beads);

    % Store particle coordination number (# particle neighbors per particle)
    % as well as # of touching particles
    bead_neighbors       = cell(num_beads, 1);
    bead_coord_count     = zeros(num_beads, 1);
    touching_beads       = cell(num_beads, 1);
    touching_beads_count = zeros(num_beads, 1);
    for i = 1 : num_beads
        neighs_across = find(beads_adjacency(i, :) ~= 0);
        neighs_down   = find(beads_adjacency(:, i) ~= 0);
        bead_neighbors{i} = sort([neighs_across(:); neighs_down(:)]);
        bead_coord_count(i)  = numel(bead_neighbors{i});
        % check which particles are actually touching
        touch_across = find(touching_adjacency(i, :) ~= 0);
        touch_down   = find(touching_adjacency(:, i) ~= 0);
        touching_beads{i} = sort([touch_across(:); touch_down(:)]);
        touching_beads_count(i)  = numel(touching_beads{i});
    end
    bead_struct.Neighbors     = bead_neighbors;
    bead_struct.CoordCount    = bead_coord_count;
    bead_struct.TouchingBeads = touching_beads;
    bead_struct.TouchingCount = touching_beads_count;

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Classify ridges2D:');

    time_log(timeLogIdx).Name = 'Classify ridges2D';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % Create cell array to store medial axes row indices by size
    med_min     = min(medax.numBeads);
    med_max     = max(medax.numBeads);
    medax.sized = cell(med_max, 1); % store row indices
    % medax.orphans = cell(med_max, 1); % store indices not matched to ridges or peaks
    for i = med_min : med_max
        inds = find(medax.numBeads == i);
        medax.sized{i} = inds(:);
    end

    %%%%%%%**********************************************************%%%%%%%
    %%%%%*************************** LOOPS ****************************%%%%%
    %%%%%%%**********************************************************%%%%%%%

    tStart = tic;

    % `loopsBySize` maps loop length to the set of loops of that length, which is
    % stored in an N x L matrix: N loops of length L. It is important to note that
    % each row of the N x L matrix contains the unique node indices (sorted in
    % increase order) that make up each loop.
    loopsBySize = findLoops(medax, ridges2D.beadNeighbors);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'SSSR:');

    time_log(timeLogIdx).Name = 'SSSR';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    %%%%%%%**********************************************************%%%%%%%
    %%%%%************************* RIDGES 1D **************************%%%%%
    %%%%%%%**********************************************************%%%%%%%

    tStart = tic;

    % `constructRidges1D` will produce a `ridges1D` struct that contains four fields:
    %  * indices_sized = A cell array of cells, where `indices_sized{i}` contains
    %                    cells of voxels (by their linear index) that belong to a
    %                    1D-ridge that is close to `i` beads. Each cell of
    %                    `indices_sized{i}` is a unique 1D-ridge.
    %  * beads_sized   = A cell array of matrices, where `beads_sized{i}` is an
    %                    `N x i' matrix where the `r`-th row contains the `i` bead
    %                    indices that the `r`-th 1D-ridge is close to.
    %  * beads_indexed = A cell array of 1D arrays, where `beads_indexed{i}` is an
    %                    array of length `N` that contains the unique (linear) ID
    %                    of each row of `beads_sized{i}`. These IDs are commonly
    %                    referred to as the hash of the bead tuple.
    %  * num           = Total number of 1D-ridges.
    %
    % The second output of `constructRidges1D` is a logical array where the `i`-th
    % component is true if the `i`-th medax voxel wasn't matched with a 1D-ridge.
    [ridges1D, medax.nonr1Dbool] = constructRidges1D(medax, loopsBySize);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Classify ridges1D:');

    time_log(timeLogIdx).Name = 'Classify ridges1D';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    %%%%%%%*************************************%%%%%%%
    %%%%%***************** PEAKS *****************%%%%%
    %%%%%%%*************************************%%%%%%%

    % Turn off any remaining medial axes that are equidistant to 3 beads, as
    % these form edge ridges1D that are missing the 3rd ridges2D because
    % it extends outside of the void space - we ignore these ridges1D
    medax.nonr1Dbool = (medax.nonr1Dbool & medax.numBeads ~= 3);

    tStart = tic;
    % Remaining medial axes voxels should contain peaks
    peaks = struct;
    peaks.L1.indices  = medax.indices(medax.nonr1Dbool);
    peaks.L1.beads    = medax.beads(:, medax.nonr1Dbool);
    peaks.L1.numBeads = medax.numBeads(medax.nonr1Dbool);

    % Cluster peaks that share the same neighboring beads
    % Assume peaks that are equidistant to a larger number of beads are the
    % true peak among all peaks that contain a subset of neighboring beads

    % Gather voxels based on number of equidistant beads
    min_bds = min(peaks.L1.numBeads);
    max_bds = max(peaks.L1.numBeads);
    bds_range = min_bds : max_bds;
    % peaks_sized cell array will contain cells of clustered peaks
    peaks_sized = cell(max_bds, 1);
    % for each number of equidistant beads...
    for i = bds_range
        % create array where each row refers to voxel and contains sorted beads
        i_pks1 = find(peaks.L1.numBeads == i);
        if isempty(i_pks1) continue; end
        beads_tuple = sort(getNearestBeads(peaks.L1.beads, i_pks1));
        [beads_tuple, sort_ind] = sortrowsReverse(beads_tuple');
        i_pks1 = i_pks1(sort_ind); % index into peaks.L1
        % store information
        peaks_sized{i}          = struct;
        peaks_sized{i}.indices  = cell(size(beads_tuple, 1), 1);
        peaks_sized{i}.beads    = zeros(i, size(beads_tuple, 1));
        peaks_sized{i}.pks1_ind = cell(size(beads_tuple, 1), 1);
        % view each row as a coordinate; change from tuples to indices
        mat_shape = repmat(num_beads, 1, i);
        beads_ind = all_sub2ind(mat_shape, beads_tuple); % already sorted
        % grab bins
        iter = 1;
        cell_cntr = 1;
        while iter <= length(beads_ind)
            beg_ind = iter;
            while iter < length(beads_ind) && beads_ind(iter + 1) == beads_ind(iter)
                iter = iter + 1;
            end
            end_ind  = iter;
            row_inds = i_pks1(beg_ind : end_ind); % indices into peaks.L1
            peaks_sized{i}.indices{cell_cntr}  = peaks.L1.indices(row_inds);
            peaks_sized{i}.beads(:, cell_cntr) = beads_tuple(beg_ind, :)';
            peaks_sized{i}.pks1_ind{cell_cntr} = row_inds;
            iter = iter + 1;
            cell_cntr = cell_cntr + 1;
        end
        remove_cells = cellfun(@isempty, peaks_sized{i}.indices);
        peaks_sized{i}.indices(remove_cells)  = [];
        peaks_sized{i}.beads(:, remove_cells) = [];
        peaks_sized{i}.pks1_ind(remove_cells) = [];
    end

    % check if smaller set of equidistant beads are contained within larger
    % set. If so, move all voxels from smaller set to ridges1D
    % Not all peaks_sized will contain voxels, so need to check for this
    nonempty_clust = find(~cellfun(@isempty, peaks_sized));
    nonempty_clust = nonempty_clust(~arrayfun(@(x) isempty(peaks_sized{x}.indices), nonempty_clust));
    removed = [];
    if length(nonempty_clust) > 1
        for h = 1 : length(nonempty_clust) - 1
            i = nonempty_clust(h);
            mat_shape = repmat(num_beads, 1, i);
            % Convert cluster i beads from tuples to index
            i_subs_ind = all_sub2ind(mat_shape, peaks_sized{i}.beads'); % already sorted
            % keep track of which rows are removed
            remove_these = [];
            for hh = h + 1 : length(nonempty_clust)
                j = nonempty_clust(hh);
                % Get all permutations of smaller set size using larger set to
                % find which smaller sets are contained in larger set of beads
                for k = 1 : size(peaks_sized{j}.beads, 2)
                    % get all permutations
                    subsets = sortrowsReverse(nchoosek(peaks_sized{j}.beads(:, k), i));
                    % change from tuples to indices
                    subsets_ind = all_sub2ind(mat_shape, subsets); % already sorted
                    % find if any i_subs_ind match subsets_ind
                    valid_rows = fastIntersect(i_subs_ind, subsets_ind, 'indices');

                    if ~isempty(valid_rows)
                        % determine how many beads to check
                        if i - 1 <= length(ridges1D.beads_indexed) && ...
                           ~isempty(ridges1D.beads_indexed{i - 1})
                            this_many = i - 1;
                        elseif i <= length(ridges1D.beads_indexed)
                            this_many = i;
                        else
                            this_many = [];
                            ii = i - 2;
                            while isempty(this_many)
                                if ~isempty(ridges1D.beads_indexed{ii})
                                    this_many = ii;
                                else
                                    ii = ii - 1;
                                end
                                if ii == 0
                                    error('ridges1D are not associated with any beads.')
                                end
                            end
                        end
                        mat_shape = repmat(num_beads, 1, this_many);
                        % find row of index in peaks.L1.indices to grab nearest
                        % beads to match to ridges1D
                        these_clusters  = peaks_sized{i}.pks1_ind(valid_rows);
                        these_pks1_ind  = cell2mat(these_clusters);
                        % associate each voxel to ridge1D and add to ridge
                        for m = 1 : length(these_pks1_ind)
                            % keep expanding nearest beads
                            numbeads2include = this_many;
                            this_ridge = -1;

                            thisPeakIdx = these_pks1_ind(m);
                            thisVoxelIdx = peaks.L1.indices(thisPeakIdx);
                            allNearBeads = getNearestBeads(peaks.L1.beads, thisPeakIdx);
                            allDists = zeros(size(allNearBeads));
                            for ii = 1 : numel(allDists)
                                allDists(ii) = sparseBeadEDTs{allNearBeads(ii)}(thisVoxelIdx, 1);
                            end
                            % Sorting distances with bead indices as well to help break
                            % ties in a deterministic way, i.e., if two beads are tied
                            % for distance, the bead with the smaller index wins.
                            [~, sortPerm] = sortrows([allDists, allNearBeads]);
                            allNearBeads = allNearBeads(sortPerm);

                            while this_ridge < 0 && numbeads2include <= peaks.L1.numBeads(these_pks1_ind(m))
                                % get all permutations
                                beadIdxs = sort(allNearBeads(1 : numbeads2include));
                                perms = sortrowsReverse(nchoosek(beadIdxs, this_many));
                                % change from tuples to indices
                                perms_ind = all_sub2ind(mat_shape, perms); % already sorted
                                % find ridge match to any permutation
                                found_match = false;
                                r_cnt = 1;
                                while ~found_match && r_cnt <= length(ridges1D.beads_indexed{this_many})
                                    if binarySearch(perms_ind, ridges1D.beads_indexed{this_many}(r_cnt)) > 0
                                        ridges1D.indices_sized{this_many}{r_cnt} = sort([ridges1D.indices_sized{this_many}{r_cnt}; ...
                                                                                         peaks.L1.indices(these_pks1_ind(m))]);
                                        this_ridge  = r_cnt;
                                        found_match = true;
                                    end
                                    r_cnt = r_cnt + 1;
                                end
                                numbeads2include = numbeads2include + 1;
                            end
                            % if no match exists, voxel is not stored
                        end
                        mat_shape = repmat(num_beads, 1, i);
                    end
                    % remove these rows from i_subs_ind
                    i_subs_ind(valid_rows) = -1;
                    % keep track of which rows to remove
                    remove_these = [remove_these; valid_rows];
                end
            end
            remove_these = unique(remove_these);
            removed = [removed; peaks_sized{i}.indices(remove_these)];
            peaks_sized{i}.indices(remove_these)  = [];
            peaks_sized{i}.beads(:, remove_these) = [];
            peaks_sized{i}.pks1_ind(remove_these) = [];
        end
    end
    % combine peak clusters
    nonempty_pks = find(~cellfun(@isempty, peaks_sized));
    tot_clusts = sum(arrayfun(@(x) length(peaks_sized{x}.indices), nonempty_pks));
    peaks.L1.clusters = cell(tot_clusts, 1);
    cntr = 1;
    for i = nonempty_pks(end) : -1 : nonempty_pks(1)
        if isempty(peaks_sized{i}), continue; end
        for j = 1 : length(peaks_sized{i}.indices)
            peaks.L1.clusters{cntr} = [peaks_sized{i}.indices{j}, ...
                                       peaks_sized{i}.pks1_ind{j}];
            cntr = cntr + 1;
        end
    end

    % For each cluster, find the peak with the max edt
    % Initialize new set of peaks, as well as updated numBeads vec
    peaks.L2.indices  = zeros(length(peaks.L1.clusters), 1);
    peaks.L2.beads    = cell(length(peaks.L1.clusters), 1);
    peaks.L2.numBeads = zeros(length(peaks.L1.clusters), 1);
    for i = 1 : length(peaks.L1.clusters)
        % if cluster is a single peak, just select it
        if size(peaks.L1.clusters{i}, 1) == 1
            peaks.L2.indices(i) = peaks.L1.clusters{i}(1, 1);
            % j is index in peaks.L1
            j = peaks.L1.clusters{i}(1, 2);
            peaks.L2.numBeads(i) = peaks.L1.numBeads(j);
            current_beads = sort(getNearestBeads(peaks.L1.beads, j));
            peaks.L2.beads{i} = current_beads(:);
        else
            % find max among peaks
            [~, max_ind] = max(data.EDT(peaks.L1.clusters{i}(:, 1)));
            peaks.L2.indices(i) = peaks.L1.clusters{i}(max_ind, 1);
            % j is index in peaks.L1
            j = peaks.L1.clusters{i}(max_ind, 2);
            peaks.L2.numBeads(i) = peaks.L1.numBeads(j);
            current_beads = sort(getNearestBeads(peaks.L1.beads, j));
            peaks.L2.beads{i} = current_beads(:);

            % Move non-peaks to appropriate ridges1D
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = [1 : max_ind - 1, max_ind + 1 : size(peaks.L1.clusters{i}, 1)]
                l = peaks.L1.clusters{i}(k, 2);
                pks_numbds = peaks.L1.numBeads(l);
                % determine how many beads to check
                if pks_numbds - 1 <= length(ridges1D.beads_indexed) && ...
                   ~isempty(ridges1D.beads_indexed{pks_numbds - 1})
                    this_many = pks_numbds - 1;
                elseif pks_numbds <= length(ridges1D.beads_indexed)
                    this_many = pks_numbds;
                else
                    this_many = [];
                    ii = pks_numbds - 2;
                    while isempty(this_many)
                        if ~isempty(ridges1D.beads_indexed{ii})
                            this_many = ii;
                        else
                            ii = ii - 1;
                        end
                        if ii == 0
                            error('ridges1D are not associated with any beads.')
                        end
                    end
                end
                mat_shape = repmat(num_beads, 1, this_many);

                % expand nearest beads looking for match
                numbeads2include = this_many;
                this_ridge = -1;

                thisVoxelIdx = peaks.L1.indices(l);
                allNearBeads = getNearestBeads(peaks.L1.beads, l);
                allDists = zeros(size(allNearBeads));
                for ii = 1 : numel(allDists)
                    allDists(ii) = sparseBeadEDTs{allNearBeads(ii)}(thisVoxelIdx, 1);
                end
                % Sorting distances with bead indices as well to help break
                % ties in a deterministic way, i.e., if two beads are tied
                % for distance, the bead with the smaller index wins.
                [~, sortPerm] = sortrows([allDists, allNearBeads]);
                allNearBeads = allNearBeads(sortPerm);

                while this_ridge < 0 && numbeads2include <= peaks.L1.numBeads(l)
                    % get all permutations
                    beadIdxs = sort(allNearBeads(1 : numbeads2include));
                    perms = sortrowsReverse(nchoosek(beadIdxs, this_many));
                    % change from tuples to indices
                    perms_ind = all_sub2ind(mat_shape, perms); % already sorted
                    % find ridge match to any permutation
                    found_match = false;
                    r_cnt = 1;
                    while ~found_match && r_cnt <= length(ridges1D.beads_indexed{this_many})
                        if binarySearch(perms_ind, ridges1D.beads_indexed{this_many}(r_cnt)) > 0
                            ridges1D.indices_sized{this_many}{r_cnt} = sort([ridges1D.indices_sized{this_many}{r_cnt}; ...
                                                                             peaks.L1.clusters{i}(k, 1)]);
                            this_ridge  = r_cnt;
                            found_match = true;
                        end
                        r_cnt = r_cnt + 1;
                    end
                    numbeads2include = numbeads2include + 1;
                end
                % if no match exists, voxel is removed
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end

    [peaks.L2.indices, order] = sort(peaks.L2.indices);
    peaks.L2.beads            = peaks.L2.beads(order);
    peaks.L2.numBeads         = peaks.L2.numBeads(order);
    peaks.L2.num              = length(peaks.L2.indices);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Peaks clustering:');

    time_log(timeLogIdx).Name = 'Peaks clustering';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    tStart = tic;
    % Cluster peaks that are PHYSICALLY TOUCHING
    touching_peaks = zeros(shape);
    touching_peaks(peaks.L2.indices) = 1;
    peaks.L2.clusters = bwconncomp(touching_peaks);
    numClust = length(peaks.L2.clusters.PixelIdxList);
    peaks.L3.indices  = zeros(numClust, 1);
    peaks.L3.beads    = cell(numClust, 1);
    peaks.L3.numBeads = zeros(numClust, 1);

    % find peak among clusters
    for i = 1 : numClust
        if length(peaks.L2.clusters.PixelIdxList{i}) > 1
            these_pks2_ind = fastIntersect(peaks.L2.indices, peaks.L2.clusters.PixelIdxList{i}, 'indices');
            [~, max_ind] = max(data.EDT(peaks.L2.clusters.PixelIdxList{i}));
            peaks.L3.indices(i)  = peaks.L2.clusters.PixelIdxList{i}(max_ind);
            peaks.pr1D_key(i, 1) = peaks.L3.indices(i);
            % collect all beads associated with all peaks in cluster
            peaks.L3.beads{i}    = unique(cell2mat(peaks.L2.beads(these_pks2_ind)));
            peaks.L3.numBeads(i) = length(peaks.L3.beads{i});

            % Move non-peaks to appropriate ridges1D
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = [1 : max_ind - 1, max_ind + 1 : length(peaks.L2.clusters.PixelIdxList{i})]
                pks_numbds = peaks.L2.numBeads(these_pks2_ind(k));
                l = binarySearch(peaks.L1.indices, peaks.L2.indices(these_pks2_ind(k)));
                % determine how many beads to check
                if pks_numbds - 1 <= length(ridges1D.beads_indexed) && ...
                   ~isempty(ridges1D.beads_indexed{pks_numbds - 1})
                    this_many = pks_numbds - 1;
                elseif pks_numbds <= length(ridges1D.beads_indexed)
                    this_many = pks_numbds;
                else
                    this_many = [];
                    ii = pks_numbds - 2;
                    while isempty(this_many)
                        if ~isempty(ridges1D.beads_indexed{ii})
                            this_many = ii;
                        else
                            ii = ii - 1;
                        end
                    end
                end
                mat_shape = repmat(num_beads, 1, this_many);

                % keep expanding nearest beads until you find a match
                numbeads2include = this_many;
                this_ridge = -1;

                thisVoxelIdx = peaks.L1.indices(l);
                allNearBeads = getNearestBeads(peaks.L1.beads, l);
                allDists = zeros(size(allNearBeads));
                for ii = 1 : numel(allDists)
                    allDists(ii) = sparseBeadEDTs{allNearBeads(ii)}(thisVoxelIdx, 1);
                end
                % Sorting distances with bead indices as well to help break
                % ties in a deterministic way, i.e., if two beads are tied
                % for distance, the bead with the smaller index wins.
                [~, sortPerm] = sortrows([allDists, allNearBeads]);
                allNearBeads = allNearBeads(sortPerm);

                while this_ridge < 0 && numbeads2include <= peaks.L1.numBeads(l)
                    % get all permutations
                    beadIdxs = sort(allNearBeads(1 : numbeads2include));
                    perms = sortrowsReverse(nchoosek(beadIdxs, this_many));
                    % change from tuples to indices
                    perms_ind = all_sub2ind(mat_shape, perms); % already sorted
                    % find ridge match to any permutation
                    found_match = false;
                    r_cnt = 1;
                    while ~found_match && r_cnt <= length(ridges1D.beads_indexed{this_many})
                        if binarySearch(perms_ind, ridges1D.beads_indexed{this_many}(r_cnt)) > 0
                            ridges1D.indices_sized{this_many}{r_cnt} = sort([ridges1D.indices_sized{this_many}{r_cnt}; ...
                                                                             peaks.L2.indices(these_pks2_ind(k))]);
                            this_ridge  = r_cnt;
                            found_match = true;
                        end
                        r_cnt = r_cnt + 1;
                    end
                    numbeads2include = numbeads2include + 1;
                end
                % if no match exists, voxel is removed
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            pk_row = binarySearch(peaks.L2.indices, peaks.L2.clusters.PixelIdxList{i});
            peaks.L3.indices(i)  = peaks.L2.clusters.PixelIdxList{i};
            peaks.pr1D_key(i, 1) = peaks.L3.indices(i);
            peaks.L3.beads{i}    = peaks.L2.beads{pk_row};
            peaks.L3.numBeads(i) = peaks.L2.numBeads(pk_row);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % with the additional beads added to updated peaks, re-check whether any
    % other peaks' beads are subsets of updated beads
    min_bds = min(peaks.L3.numBeads);
    max_bds = max(peaks.L3.numBeads);
    bds_rangelog = arrayfun(@(x) ~isempty(find(peaks.L3.numBeads == x, 1)), ...
                            min_bds : max_bds);
    bds_range = min_bds : max_bds;
    bds_range = bds_range(bds_rangelog);

    % peaks_sized cell array will contain cells of clustered peaks
    peaks_sized = cell(max_bds, 1);
    % for each number of equidistant beads...
    for i = bds_range
        % create array where each row refers to voxel and contains sorted beads
        i_pks3 = find(peaks.L3.numBeads == i);
        beads_tuple = cell2mat(peaks.L3.beads(i_pks3)')';
        [beads_tuple, sort_ind] = sortrowsReverse(beads_tuple);
        i_pks3 = i_pks3(sort_ind); % index into peaks.L3
        % store information
        peaks_sized{i}          = struct;
        peaks_sized{i}.indices  = cell(size(beads_tuple, 1), 1); % idk why this is a cell...
        peaks_sized{i}.beads    = zeros(i, size(beads_tuple, 1));
        peaks_sized{i}.pks3_ind = cell(size(beads_tuple, 1), 1);
        % view each row as a coordinate; change from tuples to indices
        mat_shape = repmat(num_beads, 1, i);
        beads_ind = all_sub2ind(mat_shape, beads_tuple); % already sorted
        % grab bins
        iter = 1;
        cell_cntr = 1;
        while iter <= length(beads_ind)
            beg_ind = iter;
            while iter < length(beads_ind) && beads_ind(iter + 1) == beads_ind(iter)
                iter = iter + 1;
            end
            end_ind  = iter;
            row_inds = i_pks3(beg_ind : end_ind); % indices into peaks.L3
            peaks_sized{i}.indices{cell_cntr}  = peaks.L3.indices(row_inds);
            peaks_sized{i}.beads(:, cell_cntr) = beads_tuple(beg_ind, :)';
            peaks_sized{i}.pks3_ind{cell_cntr} = row_inds;
            iter = iter + 1;
            cell_cntr = cell_cntr + 1;
        end
        remove_cells = cellfun(@isempty, peaks_sized{i}.indices);
        peaks_sized{i}.indices(remove_cells)  = [];
        peaks_sized{i}.beads(:, remove_cells) = [];
        peaks_sized{i}.pks3_ind(remove_cells) = [];
    end

    % check if smaller set of equidistant beads are contained within larger
    % set. If so, move all voxels from smaller set to ridges1D
    nonempty_clust = find(~cellfun(@isempty, peaks_sized));
    nonempty_clust = nonempty_clust(~arrayfun(@(x) isempty(peaks_sized{x}.indices), nonempty_clust));
    check_these = [];
    removed = [];
    if length(nonempty_clust) > 1
        for h = 1 : length(nonempty_clust) - 1
            i = nonempty_clust(h);
            mat_shape = repmat(num_beads, 1, i);
            % Convert cluster i beads from tuples to index
            i_subs_ind = all_sub2ind(mat_shape, peaks_sized{i}.beads'); % already sorted
            % keep track of which rows are removed
            remove_these = [];
            for hh = h + 1 : length(nonempty_clust)
                j = nonempty_clust(hh);
                % Get all permutations of smaller set size using larger set to
                % find which smaller sets are contained in larger set of beads
                for k = 1 : size(peaks_sized{j}.beads, 2)
                    % get all permutations
                    subsets = sortrowsReverse(nchoosek(peaks_sized{j}.beads(:, k), i));
                    % change from tuples to indices
                    subsets_ind = all_sub2ind(mat_shape, subsets); % already sorted
                    % find if any i_subs_ind match subsets_ind
                    % For i_subs_row, each row refers to bead cluster; values are
                    % irrelevant since they all refer to the peaks_sized{j}.beads(k) cluster
                    i_subs_row = arrayfun(@(x) binarySearch(subsets_ind, x), i_subs_ind);
                    valid_rows = find(i_subs_row >= 1);

                    if ~isempty(valid_rows)
                        % determine how many beads to check for ridge
                        if i - 1 <= length(ridges1D.beads_indexed) && ...
                           ~isempty(ridges1D.beads_indexed{i - 1})
                            this_many = i - 1;
                        elseif i <= length(ridges1D.beads_indexed)
                            this_many = i;
                        else
                            this_many = [];
                            ii = i - 2;
                            while isempty(this_many)
                                if ~isempty(ridges1D.beads_indexed{ii})
                                    this_many = ii;
                                else
                                    ii = ii - 1;
                                end
                            end
                        end
                        mat_shape = repmat(num_beads, 1, this_many);
                        % find row of index in peaks.L1.indices to grab nearest
                        % beads to match to ridges1D
                        these_clusters = peaks_sized{i}.pks3_ind(valid_rows);
                        these_inds     = peaks.L3.indices(cell2mat(these_clusters));
                        these_pks1_ind = fastIntersect(peaks.L1.indices, sort(these_inds), 'indices');
                        % Associate each voxel to ridge1D and add to ridge
                        for m = 1 : length(these_pks1_ind)
                            % Keep expanding nearest beads until you find a match
                            numbeads2include = this_many;
                            this_ridge = -1;

                            thisPeakIdx = these_pks1_ind(m);
                            thisVoxelIdx = peaks.L1.indices(thisPeakIdx);
                            allNearBeads = getNearestBeads(peaks.L1.beads, thisPeakIdx);
                            allDists = zeros(size(allNearBeads));
                            for ii = 1 : numel(allDists)
                                allDists(ii) = sparseBeadEDTs{allNearBeads(ii)}(thisVoxelIdx, 1);
                            end
                            % Sorting distances with bead indices as well to help break
                            % ties in a deterministic way, i.e., if two beads are tied
                            % for distance, the bead with the smaller index wins.
                            [~, sortPerm] = sortrows([allDists, allNearBeads]);
                            allNearBeads = allNearBeads(sortPerm);

                            while this_ridge < 0 && numbeads2include <= peaks.L1.numBeads(these_pks1_ind(m))
                                % get all permutations
                                beadIdxs = sort(allNearBeads(1 : numbeads2include));
                                perms = sortrowsReverse(nchoosek(beadIdxs, this_many));
                                % change from tuples to indices
                                perms_ind = all_sub2ind(mat_shape, perms); % already sorted
                                % find ridge match to any permutation
                                found_match = false;
                                r_cnt = 1;
                                while ~found_match && r_cnt <= length(ridges1D.beads_indexed{this_many})
                                    if binarySearch(perms_ind, ridges1D.beads_indexed{this_many}(r_cnt)) > 0
                                        ridges1D.indices_sized{this_many}{r_cnt} = sort([ridges1D.indices_sized{this_many}{r_cnt}; ...
                                                                                         peaks.L1.indices(these_pks1_ind(m))]);
                                        this_ridge  = r_cnt;
                                        found_match = true;
                                    end
                                    r_cnt = r_cnt + 1;
                                end
                                numbeads2include = numbeads2include + 1;
                            end
                            % if no match exists, medial axis voxel is removed
                            check_these = [check_these; peaks.L1.indices(these_pks1_ind(m))];
                        end
                        mat_shape = repmat(num_beads, 1, i);
                    end
                    % remove these rows from i_subs_ind
                    i_subs_ind(valid_rows) = -1;
                    % keep track of which rows to remove
                    remove_these = [remove_these; valid_rows];
                end
            end
            remove_these = unique(remove_these);
            removed = [removed; peaks_sized{i}.indices(remove_these)];
            peaks_sized{i}.indices(remove_these)  = [];
            peaks_sized{i}.beads(:, remove_these) = [];
            peaks_sized{i}.pks3_ind(remove_these) = [];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % combine peak clusters
    nonempty_pks = find(~cellfun(@isempty, peaks_sized));
    tot_clusts = sum(arrayfun(@(x) length(peaks_sized{x}.indices), nonempty_pks));
    peaks.L3.clusters = cell(tot_clusts, 1);
    cntr = 1;
    for i = nonempty_pks(end) : -1 : nonempty_pks(1)
        if ~isempty(peaks_sized{i})
            for j = 1 : length(peaks_sized{i}.indices)
                peaks.L3.clusters{cntr} = [peaks_sized{i}.indices{j}, ...
                                           peaks_sized{i}.pks3_ind{j}];
                cntr = cntr + 1;
            end
        end
    end
    % For each cluster, find the peak with the max edt
    % Initialize new set of peaks, as well as updated numBeads vec
    peaks.L4.indices  = zeros(length(peaks.L3.clusters), 1);
    peaks.L4.beads    = cell(length(peaks.L3.clusters), 1);
    peaks.L4.numBeads = zeros(length(peaks.L3.clusters), 1);
    for i = 1 : length(peaks.L3.clusters)
        % if cluster is a single peak, just select it
        if size(peaks.L3.clusters{i}, 1) == 1
            peaks.L4.indices(i) = peaks.L3.clusters{i}(1, 1);
            % j is index in peaks.L3
            j = peaks.L3.clusters{i}(1, 2);
            peaks.L4.numBeads(i) = peaks.L3.numBeads(j);
            peaks.L4.beads{i}    = peaks.L3.beads{j};
        else
            % find max among peaks
            [~, max_ind] = max(data.EDT(peaks.L3.clusters{i}(:, 1)));
            peaks.L4.indices(i) = peaks.L3.clusters{i}(max_ind, 1);
            % j is index in peaks.L1
            j = peaks.L3.clusters{i}(max_ind, 2);
            peaks.L4.numBeads(i) = peaks.L3.numBeads(j);
            peaks.L4.beads{i}    = peaks.L3.beads{j};

            % Move non-peaks to appropriate ridges1D
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = [1 : max_ind - 1, max_ind + 1 : size(peaks.L3.clusters{i}, 1)]
                l = binarySearch(peaks.L1.indices, peaks.L3.clusters{i}(k, 1));
                pks_numbds = peaks.L3.numBeads(peaks.L3.clusters{i}(k, 2));
                % determine how many beads to check
                if pks_numbds - 1 <= length(ridges1D.beads_indexed) && ...
                   ~isempty(ridges1D.beads_indexed{pks_numbds - 1})
                    this_many = pks_numbds - 1;
                elseif pks_numbds <= length(ridges1D.beads_indexed)
                    this_many = pks_numbds;
                else
                    this_many = [];
                    ii = pks_numbds - 2;
                    while isempty(this_many)
                        if ~isempty(ridges1D.beads_indexed{ii})
                            this_many = ii;
                        else
                            ii = ii - 1;
                        end
                    end
                end
                mat_shape = repmat(num_beads, 1, this_many);

                % keep expanding nearest beads until you find a match
                numbeads2include = this_many;
                this_ridge = -1;

                thisVoxelIdx = peaks.L1.indices(l);
                allNearBeads = getNearestBeads(peaks.L1.beads, l);
                allDists = zeros(size(allNearBeads));
                for ii = 1 : numel(allDists)
                    allDists(ii) = sparseBeadEDTs{allNearBeads(ii)}(thisVoxelIdx, 1);
                end
                % Sorting distances with bead indices as well to help break
                % ties in a deterministic way, i.e., if two beads are tied
                % for distance, the bead with the smaller index wins.
                [~, sortPerm] = sortrows([allDists, allNearBeads]);
                allNearBeads = allNearBeads(sortPerm);

                while this_ridge < 0 && numbeads2include <= peaks.L1.numBeads(l)
                    % get all permutations
                    beadIdxs = sort(allNearBeads(1 : numbeads2include));
                    perms = sortrowsReverse(nchoosek(beadIdxs, this_many));
                    % change from tuples to indices
                    perms_ind = all_sub2ind(mat_shape, perms); % already sorted
                    % find ridge match to any permutation
                    found_match = false;
                    r_cnt = 1;
                    while ~found_match && r_cnt <= length(ridges1D.beads_indexed{this_many})
                        if binarySearch(perms_ind, ridges1D.beads_indexed{this_many}(r_cnt)) > 0
                            ridges1D.indices_sized{this_many}{r_cnt} = sort([ridges1D.indices_sized{this_many}{r_cnt}; ...
                                                                             peaks.L3.clusters{i}(k, 1)]);
                            this_ridge  = r_cnt;
                            found_match = true;
                        end
                        r_cnt = r_cnt + 1;
                    end
                    numbeads2include = numbeads2include + 1;
                end
                % if no match exists, voxel is removed
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end

    [peaks.L4.indices, order] = sort(peaks.L4.indices);
    peaks.L4.beads            = peaks.L4.beads(order);
    peaks.L4.numBeads         = peaks.L4.numBeads(order);
    peaks.L4.num              = length(peaks.L4.indices);
    % Store EDT of peaks
    peaks.L4.EDT = data.EDT(peaks.L4.indices);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Peaks touching:');

    time_log(timeLogIdx).Name = 'Peaks touching';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % combine ridges1D
    % store indices and beads
    ridges1D.indices = cell(ridges1D.num, 1);
    ridges1D.beads   = cell(ridges1D.num, 1);
    % create matrix look-up key
    max_length_j     = max(cellfun(@length, ridges1D.indices_sized));
    r1D_sizedkey     = zeros(ridges1D.num, 2);
    cntr = 1;
    for i = 1 : length(ridges1D.indices_sized)
        for j = 1 : length(ridges1D.indices_sized{i})
            ridges1D.indices{cntr} = ridges1D.indices_sized{i}{j};
            ridges1D.beads{cntr}   = ridges1D.beads_sized{i}(j, :);
            r1D_sizedkey(cntr, :)  = [i, j];
            cntr = cntr + 1;
        end
    end
    % create sparse matrix for look-up key
    r1D_sizedkey = sparse(r1D_sizedkey(:, 1), r1D_sizedkey(:, 2), ...
                          1 : size(r1D_sizedkey, 1), ...
                          length(ridges1D.indices_sized), max_length_j);

    tStart = tic;
    % Associate peaks and ridges
    % peaks.pr1D_key:
    %       Column 1:    peak indices
    %       Column 2-16: associated ridge #s
    % ridges1D.rp_key:
    %       Rows:        each row refers to a unique ridge1D (in order)
    %       Columns:     associated peaks
    peaks.pr1D_key       = zeros(length(peaks.L4.indices), 16);
    peaks.pr1D_key(:, 1) = peaks.L4.indices;
    peaks.pr1D_counter   = ones(size(peaks.pr1D_key, 1), 1);

    ridges1D.rp_key      = zeros(ridges1D.num, 10);
    ridges1D.rp_counter  = zeros(ridges1D.num, 1);

    for i = 1 : peaks.L4.num
        pks_numbds = peaks.L4.numBeads(i);
        for j = pks_numbds : -1 : 3
            % determine how many beads to check
            if j <= length(ridges1D.beads_indexed) && ~isempty(ridges1D.beads_indexed{j})
                mat_shape = repmat(num_beads, 1, j);

                % get all permutations
                perms = sortrowsReverse(nchoosek(peaks.L4.beads{i}, j));
                % change from tuples to indices
                perms_ind = all_sub2ind(mat_shape, perms); % already sorted
                % find ridge matches to any permutation
                for k = 1 : length(perms_ind)
                    this_ridge_ind = binarySearch(ridges1D.beads_indexed{j}, perms_ind(k));
                    if this_ridge_ind > 0
                        this_ridge = r1D_sizedkey(j, this_ridge_ind);
                        % check to make sure ridge isn't already
                        % associated with peak
                        if ~fastIntersect(ridges1D.rp_key(this_ridge, 1 : ridges1D.rp_counter(this_ridge)), ...
                                          peaks.L4.indices(i), 'bool')
                            ridges1D.rp_counter(this_ridge) = ridges1D.rp_counter(this_ridge) + 1;
                            ridges1D.rp_key(this_ridge, ridges1D.rp_counter(this_ridge)) = peaks.L4.indices(i);
                            % keep rows sorted
                            ridges1D.rp_key(this_ridge, 1 : ridges1D.rp_counter(this_ridge)) = ...
                                            sort(ridges1D.rp_key(this_ridge, 1 : ridges1D.rp_counter(this_ridge)));

                            peaks.pr1D_counter(i) = peaks.pr1D_counter(i) + 1;
                            peaks.pr1D_key(i, peaks.pr1D_counter(i)) = this_ridge;
                            peaks.pr1D_key(i, 2 : peaks.pr1D_counter(i)) = sort(peaks.pr1D_key(i, 2 : peaks.pr1D_counter(i)));
                        end
                    end
                end
            end
        end
    end
    peaks.pr1D_counter = peaks.pr1D_counter - 1;

    % Store different ridge1D types
    ridges1D.pks0 = find(ridges1D.rp_counter == 0);
    ridges1D.pks1 = find(ridges1D.rp_counter == 1);
    ridges1D.pks2 = find(ridges1D.rp_counter >= 2);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Filling peaks ridges keys:');

    time_log(timeLogIdx).Name = 'Filling peaks ridges keys';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    %%%%%%%*************************************%%%%%%%
    %%%%%***************** DOORS *****************%%%%%
    %%%%%%%*************************************%%%%%%%

    tStart = tic;
    % Locate doors by finding local mins of unique ridges
    % average location of lowest voxels in ridge1D
    ridges1D.doors      = cell(ridges1D.num, 1);
    ridges1D.lowest_avg = zeros(ridges1D.num, 3);
    % index of voxel lying closest to average location
    ridges1D.min = zeros(ridges1D.num, 1);
    for i = 1 : ridges1D.num
        % many ridges1Ds have multiple 'lowest' EDT voxels
        lowest_voxs = data.EDT(ridges1D.indices{i}) == min(data.EDT(ridges1D.indices{i}));
        ridges1D.lowest_avg(i, :) = mean(voxels(ridges1D.indices{i}(lowest_voxs), :), 1);
        % use the voxel closest to the average distance of all lowest voxels
        avg_idx = knnsearch(voxels(ridges1D.indices{i}, :), ridges1D.lowest_avg(i, :));
        ridges1D.min(i) = ridges1D.indices{i}(avg_idx);
    end

    % per min voxel, find voxel on nearest 3 beads; each column is a
    % different bead
    ridges1D.bead_voxs = zeros(ridges1D.num, 3);
    for i = 1 : ridges1D.num
        % checking nearest 3 beads since need 3 pts to form a plane
        for h = 1 : 3
            edgebead_voxs = data.edgeIndices{ridges1D.beads{i}(h)};
            norm_vec = zeros(length(edgebead_voxs), 1);
            for j = 1 : length(edgebead_voxs)
                norm_vec(j) = norm(voxels(edgebead_voxs(j), :) - voxels(ridges1D.min(i), :));
            end
            [~, minDist_idx] = min(norm_vec);
            ridges1D.bead_voxs(i, h) = edgebead_voxs(minDist_idx);
        end
    end

    % Split ridges at door by finding plane that bisects ridge at door
    % location (avg of lowest voxels)

    % split ridge1D data is going to be stored as a cell array:
    % columns: unique ridge1Ds
    % row 1  : peak associated with upper half of ridge
    % row 2  : ridge1D indices of voxels in upper half
    % row 3  : ridge1D indices of voxels in lower half
    % row 4  : peak associated with lower half
    ridges1D.split = cell(4, ridges1D.num);
    for i = 1 : ridges1D.num
        ridges1D.doors{i} = struct;
        % find normal vector to plane
        [center, radius, v1, v2] = circlefit3d(voxels(ridges1D.bead_voxs(i, 1), :), ...
                                               voxels(ridges1D.bead_voxs(i, 2), :), ...
                                               voxels(ridges1D.bead_voxs(i, 3), :));
        if ~isempty(v1 + v2)
            nvec = cross(v1, v2);
        else
            nvec = [NaN; NaN; NaN];
        end

        % use dot product and average location of lowest EDT voxels to
        % determine whether ridge1D voxels lie above or below plane;
        % include peak voxels as well
        % ridgespeaks_ind is (x,y,z) of ridges1D stacked on (x,y,z) of peaks
        ridgepeaks_ind = [voxels(ridges1D.indices{i}, :); ...
                          voxels(ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i))', :)];
        above_below = arrayfun(@(idx) dot(nvec, ridgepeaks_ind(idx, :) - ridges1D.lowest_avg(i, :)), ...
                                      (1 : size(ridgepeaks_ind, 1)))';
        % fill ridges1D.split cell array
        rp_c = ridges1D.rp_counter(i);
        while rp_c > 0
            if above_below(end) > 0
                ridges1D.split{1, i} = [ridges1D.split{1, i}; ridges1D.rp_key(i, rp_c)];
                above_below(end) = [];
            else % lumping in < and =0
                ridges1D.split{4, i} = [ridges1D.split{4, i}; ridges1D.rp_key(i, rp_c)];
                above_below(end) = [];
            end
            rp_c = rp_c - 1;
        end
        % sort
        ridges1D.split{1, i} = sort(ridges1D.split{1, i});
        ridges1D.split{4, i} = sort(ridges1D.split{4, i});
        % split ridge1D voxels
        above = above_below > 0;
        below = above_below <= 0; % lumping in < and =0
        ridges1D.split{2, i} = sort(ridges1D.indices{i}(above));
        ridges1D.split{3, i} = sort(ridges1D.indices{i}(below));

        % store door data per ridge
        ridges1D.doors{i}.center = center;
        ridges1D.doors{i}.radius = radius;
        ridges1D.doors{i}.v1     = v1;
        ridges1D.doors{i}.v2     = v2;
        % normal vector
        ridges1D.doors{i}.nvec   = nvec;
    end

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Collect Doors Data:');

    time_log(timeLogIdx).Name = 'Collect doors data';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    %%%%%%%**************************************************%%%%%%%
    %%%%%***************** ADD & REMOVE PEAKS *****************%%%%%
    %%%%%%%**************************************************%%%%%%%

    % To ensure an accurate medial axes network, we force all interior ridges1D to
    % be flanked by 2 peaks. This may require us to add or remove peaks as
    % needed.

    tStart = tic;
    % Find edge ridges1D
    ridges1D.edge_bool = false(ridges1D.num, 1);
    % check which ridges1D flanked by < 2 peaks intersect the edge
    thicker_edge = edgeVoidVoxels(bead_voxels, vsVoxels, shape, dx, 2);
    for i = [ridges1D.pks0; ridges1D.pks1]'
        if fastIntersect(thicker_edge, ridges1D.indices{i}, 'bool')
            ridges1D.edge_bool(i) = true;
        end
    end

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Find edge ridges:');

    time_log(timeLogIdx).Name = 'Find edge ridges';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    tStart = tic;
    % MOVING PEAKS IN SPLIT
    % As a first step: for any ridge flanked by 2 peaks on one end and 0
    % peaks on the other, split peaks to flank both ends of ridge
    one_end_bool = false(ridges1D.num, 1);
    for i = 1 : ridges1D.num
        if ~ridges1D.edge_bool(i) && (isempty(ridges1D.split{1, i}) || isempty(ridges1D.split{4, i}))
            if length(ridges1D.split{1, i}) > 1 || length(ridges1D.split{4, i}) > 1
                one_end_bool(i) = true;
            end
        end
    end
    split_these = find(one_end_bool);
    for i = split_these(:)'
        if ~isempty(ridges1D.split{1, i})
            this_end  = 1;
            other_end = 4;
        elseif ~isempty(ridges1D.split{4, i})
            this_end  = 4;
            other_end = 1;
        end
        if length(ridges1D.split{this_end, i}) > 2
            pairs = nchoosek(ridges1D.split{this_end, i}, 2);
            % choose peak that is farthest from other peaks
            list_dist = arrayfun(@(x) vecnorm(voxels(pairs(x, 1), :) - ...
                                              voxels(pairs(x, 2), :), 2, 2), ...
                                 1 : size(pairs, 1));
            % overkill in case you have more than 3 pks
            keep_bool = false(length(ridges1D.split{this_end, i}), 1);
            move_ind = 1 : length(ridges1D.split{this_end, i});
            tracker = 0;
            while length(move_ind) > 1
                % keep the peaks closest to each other together; weed out
                % the farthest one
                [~, keep_tog] = min(list_dist);
                keep_bool(fastIntersect(ridges1D.split{this_end, i}, pairs(keep_tog, :), ...
                          'bool vec')) = true;
                move_ind = find(~keep_bool);
                if length(move_ind) > 1
                    list_dist(keep_tog) = [];
                    pairs(keep_tog, :)  = [];
                elseif tracker > 100
                    error('Something wrong with the logic here')
                end
                tracker = tracker + 1;
            end
            % move peak
            ridges1D.split{other_end, i} = ridges1D.split{this_end, i}(move_ind);
            ridges1D.split{this_end, i}(move_ind) = [];
        else
            % move peak that is closest to normal plane
            dist2plane = arrayfun(@(x) norm(voxels(ridges1D.split{this_end, i}(x), :) - ...
                                            ridges1D.lowest_avg(i, :)), ...
                                  1 : length(ridges1D.split{this_end, i}));
            [~, min_ind] = min(dist2plane);
            ridges1D.split{other_end, i} = ridges1D.split{this_end, i}(min_ind);
            ridges1D.split{this_end, i}(min_ind) = [];
        end
    end
    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Moving peaks flanking:');

    time_log(timeLogIdx).Name = 'Moving peaks flanking';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % REMOVING PEAKS
    % Ensure no ridge1D is flanked by more than 2 peaks
    tStart = tic;
    ridges1D_removepks = find(ridges1D.rp_counter >= 3);
    r1D_2pks = find(ridges1D.rp_counter == 2);

    % Next peaks level:
    peaks.L5 = peaks.L4;

    for i = ridges1D_removepks(:)'
        % check both ends for multiple peaks
        for j = [1, 4]
            check_pks = ridges1D.split{j, i};
            if length(check_pks) > 1 && length(check_pks) > length(unique(check_pks))
                % remove duplicates
                ridges1D.split{j, i} = unique(ridges1D.split{j, i});
                update_key = unique(ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i)));
                ridges1D.rp_counter(i) = length(update_key);
                ridges1D.rp_key(i, :) = 0;
                ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i)) = update_key;
                ridges1D.indices{i} = unique(ridges1D.indices{i});

                check_pks = ridges1D.split{j, i};
            end
            while numel(check_pks) > 1
                % check whether ridges exist between peaks. If not, combine peaks
                % using nchoosek to account for potentially > 2 pks
                pkperms  = nchoosek(check_pks, 2);
                no_match = false(size(pkperms, 1), 1);
                for k = 1 : size(pkperms, 1)
                    matches_tot = arrayfun(@(x) fastIntersect(ridges1D.rp_key(x, 1:2), ...
                                                              sort(pkperms(k, :)), 'all'), ...
                                           r1D_2pks); % only checking against ridges flanked by 2 pks
                    if sum(matches_tot) == 0
                        % confirm they are physically close, then store to combine
                        phys_dist = vecnorm(voxels(pkperms(k, 1), :) - voxels(pkperms(k, 2), :), 2, 2);
                        pk_ind1   = binarySearch(peaks.L5.indices, pkperms(k, 1));
                        pk_ind2   = binarySearch(peaks.L5.indices, pkperms(k, 2));
                        EDT_radii = peaks.L5.EDT(pk_ind1) + peaks.L5.EDT(pk_ind2);
                        if phys_dist < EDT_radii
                            no_match(k) = true;
                        end
                    end
                end
                % combine all peaks that do not share ridge1D
                combine_pks = unique(pkperms(no_match, :));
                if ~isempty(combine_pks)
                    % choose peak that is either farthest from peak on
                    % opposite end (if exists), or farthest from normal plane
                    if j == 1 && length(ridges1D.split{4, i}) == 1
                        dist_from = arrayfun(@(x) norm(voxels(combine_pks(x), :) - ...
                                                       voxels(ridges1D.split{4, i}, :)), ...
                                             1 : length(combine_pks));
                    elseif j == 4 && length(ridges1D.split{1, i}) == 1
                        dist_from = arrayfun(@(x) norm(voxels(combine_pks(x), :) - ...
                                                       voxels(ridges1D.split{1, i}, :)), ...
                                             1 : length(combine_pks));
                    else
                        dist_from = arrayfun(@(x) norm(voxels(combine_pks(x), :) - ...
                                                       ridges1D.lowest_avg(i, :)), ...
                                              1 : length(combine_pks));
                    end
                    [~, max_ind] = max(dist_from);
                    chosen_pk = combine_pks(max_ind);
                    remainder_pks = combine_pks([1 : max_ind - 1, max_ind + 1 : end]);
                    % update chosen_pks for while loop
                    rmv_inds   = fastIntersect(check_pks, remainder_pks, 'indices');
                    check_pks(rmv_inds) = [];
                    % update ridges1D key by replacing peaks with updated peak
                    for k = [1 : i - 1, i + 1 : ridges1D.num]
                        column_ind = fastIntersect(ridges1D.rp_key(k, 1 : ridges1D.rp_counter(k)), ...
                                                   remainder_pks, 'indices');
                        if ~isempty(column_ind)
                            % replace
                            if length(column_ind) > 1 && ridges1D.rp_counter(k) <= 1
                                error('A single ridge1D is flanked by 2 peaks that will be removed.')
                            else
                                % check to make sure ridge isn't already
                                % associated with chosen_pk
                                if binarySearch(ridges1D.rp_key(k, 1 : ridges1D.rp_counter(k)), ...
                                                chosen_pk) > 0
                                    % remove remainder_pks
                                    ridges1D.rp_counter(k) = ridges1D.rp_counter(k) - 1;
                                    ridges1D.rp_key(k, column_ind) = 0;
                                    ridges1D.rp_key(k, :) = shiftZeros(ridges1D.rp_key(k, :));
                                    ridges1D.rp_key(k, 1 : ridges1D.rp_counter(k)) = ...
                                        sort(ridges1D.rp_key(k, 1 : ridges1D.rp_counter(k)));

                                    if fastIntersect(ridges1D.split{1, k}, remainder_pks, 'bool')
                                        ridges1D.split{1, k}(fastIntersect(ridges1D.split{1, k}, ...
                                                             remainder_pks, 'indices')) = [];
                                    elseif fastIntersect(ridges1D.split{4, k}, remainder_pks, 'bool')
                                        ridges1D.split{4, k}(fastIntersect(ridges1D.split{4, k}, ...
                                                             remainder_pks, 'indices')) = [];
                                    else
                                        error('Ridge1D was associated with a peak in its key that does not exist in split.')
                                    end
                                else
                                    % key
                                    ridges1D.rp_key(k, column_ind) = chosen_pk;
                                    ridges1D.rp_key(k, 1 : ridges1D.rp_counter(k)) = ...
                                        sort(ridges1D.rp_key(k, 1 : ridges1D.rp_counter(k)));
                                    % split
                                    if fastIntersect(ridges1D.split{1, k}, remainder_pks, 'bool')
                                        ridges1D.split{1, k}(fastIntersect(ridges1D.split{1, k}, ...
                                                             remainder_pks, 'indices')) = chosen_pk;
                                        ridges1D.split{1, k} = sort(ridges1D.split{1, k});
                                    elseif fastIntersect(ridges1D.split{4, k}, remainder_pks, 'bool')
                                        ridges1D.split{4, k}(fastIntersect(ridges1D.split{4, k}, ...
                                                             remainder_pks, 'indices')) = chosen_pk;
                                        ridges1D.split{4, k} = sort(ridges1D.split{4, k});
                                    else
                                        error('Ridge1D was associated with a peak in its key that does not exist in split.')
                                    end
                                end
                            end
                        end
                    end
                    % remove peak from ridge i
                    these_col = fastIntersect(ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i)), ...
                                              remainder_pks, 'indices');

                    ridges1D.rp_counter(i) = ridges1D.rp_counter(i) - length(these_col);
                    ridges1D.rp_key(i, these_col) = 0;
                    ridges1D.rp_key(i, :) = shiftZeros(ridges1D.rp_key(i, :));
                    ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i)) = ...
                                        sort(ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i)));

                    if fastIntersect(ridges1D.split{1, i}, remainder_pks, 'bool')
                        ridges1D.split{1, i}(fastIntersect(ridges1D.split{1, i}, ...
                                             remainder_pks, 'indices')) = [];
                    elseif fastIntersect(ridges1D.split{4, i}, remainder_pks, 'bool')
                        ridges1D.split{4, i}(fastIntersect(ridges1D.split{4, i}, ...
                                             remainder_pks, 'indices')) = [];
                    else
                        error('Ridge1D was associated with a peak in its key that does not exist in split.')
                    end
                    % remove remainder peaks
                    if length(combine_pks) > 1
                        for k = [1 : max_ind - 1, max_ind + 1 : length(combine_pks)]
                            ind = find(peaks.L5.indices == combine_pks(k));
                            peaks.L5.indices(ind)   = [];
                            peaks.L5.beads(ind)     = [];
                            peaks.L5.numBeads(ind)  = [];
                            peaks.L5.num            = peaks.L5.num - 1;
                            peaks.L5.EDT(ind)       = [];
                            % peaks key
                            peaks.pr1D_key(ind, :)  = [];
                            peaks.pr1D_counter(ind) = [];
                        end
                        % move old peak index onto ridge
                        add2r1D = remainder_pks;
                        ridges1D.indices{i} = sort([ridges1D.indices{i}; ...
                                                    add2r1D(:)]);
                    end
                elseif length(check_pks) > 1
                    % remove extra peaks from ridge i only
                    % choose peak that is either farthest from peak on
                    % opposite end (if exists), or farthest from normal plane
                    if j == 1 && length(ridges1D.split{4, i}) == 1
                        dist_from = arrayfun(@(x) norm(voxels(check_pks(x), :) - ...
                                                       voxels(ridges1D.split{4, i}, :)), ...
                                             1 : length(check_pks));
                    elseif j == 4 && length(ridges1D.split{1, i}) == 1
                        dist_from = arrayfun(@(x) norm(voxels(check_pks(x), :) - ...
                                                       voxels(ridges1D.split{1, i}, :)), ...
                                             1 : length(check_pks));
                    else
                        dist_from = arrayfun(@(x) norm(voxels(check_pks(x), :) - ...
                                                       ridges1D.lowest_avg(i, :)), ...
                                             1 : length(check_pks));
                    end
                    [~, max_ind] = max(dist_from);
                    chosen_pk = check_pks(max_ind);
                    remainder_pks = check_pks([1 : max_ind - 1, max_ind + 1 : end]);

                    these_col = fastIntersect(ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i)), ...
                                              remainder_pks, 'indices');
                    ridges1D.rp_counter(i) = ridges1D.rp_counter(i) - length(these_col);
                    ridges1D.rp_key(i, these_col) = 0;
                    ridges1D.rp_key(i, :) = shiftZeros(ridges1D.rp_key(i, :));
                    ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i)) = ...
                        sort(ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i)));

                    these_ind = fastIntersect(ridges1D.split{j,i}, remainder_pks, 'indices');
                    ridges1D.split{j, i}(these_ind) = [];
                    % update chosen_pks for while loop
                    rmv_inds   = fastIntersect(check_pks, remainder_pks, 'indices');
                    check_pks(rmv_inds) = [];
                end
            end
        end
    end
    % Update different ridge1D types
    ridges1D.pks0 = find(ridges1D.rp_counter == 0);
    ridges1D.pks1 = find(ridges1D.rp_counter == 1);
    ridges1D.pks2 = find(ridges1D.rp_counter >= 2);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Remove peaks:');

    time_log(timeLogIdx).Name = 'Remove peaks';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % ADDING PEAKS
    % Force each interior ridge1D to be flanked by 2 peaks
    tStart = tic;
    % these are the interior ridge pts associated with less than 1 peak
    ridges1D_pks0_1  = find(ridges1D.rp_counter < 2);
    ridges1D_misspks = ridges1D_pks0_1(~ridges1D.edge_bool(ridges1D_pks0_1));

    % Next peaks level:
    peaks.L6 = peaks.L5;

    if ~isempty(ridges1D_misspks)
        % Store which ridges to remove
        remove_these_r1D = false(ridges1D.num, 1);
        for i = ridges1D_misspks(:)'
            % determine how many missing peaks
            num_miss = 2 - ridges1D.rp_counter(i);
            pks2add = zeros(num_miss, 1);
            % determine how many neighboring beads of ridge
            num_r1Dbds = length(ridges1D.beads{i});
            % first check whether all other ridges associated with existing
            % peak are flanked by 2 peaks; if so, remove ridge
            if num_miss == 1
                solo_pk = ridges1D.rp_key(i, 1);
                assoc_r1D = unique([find(ridges1D.rp_key(:,1) == solo_pk); ...
                                    find(ridges1D.rp_key(:,2) == solo_pk)]);
                assoc_r1D(assoc_r1D == i) = [];
                assoc_r1D = assoc_r1D(~ridges1D.edge_bool(assoc_r1D)); % remove edge ridges
                if fastIntersect(ridges1D.pks2, assoc_r1D, 'all')
                    % remove ridge
                    remove_these_r1D(i) = true;
                end
                            end
            % if ridge is not to be removed, proceed
            if remove_these_r1D(i) ~= true
                % collect peaks that share neighboring beads of ridge (1
                % less than ridge's number of beads)

                % The following code finds the set of all j such that any subset of size
                % num_r1Dbds - 1 of peaks.L6.beads{j} is contained in ridges1D.beads{i}.
                % However, this is equivalent to finding the set of all j such that
                % the intersection of peaks.L6.beads{j} and ridges1D.beads{i} contain
                % at least num_r1Dbds - 1 elements.

                these_pks1 = [];
                for j = 1 : peaks.L6.num
                    subset_pkbds = nchoosek(peaks.L6.beads{j}, num_r1Dbds - 1);
                    for k = 1 : size(subset_pkbds, 1)
                        if fastIntersect(ridges1D.beads{i}, subset_pkbds(k,:), 'all')
                            these_pks1 = [these_pks1; j];
                            break;
                        end
                    end
                end

                these_pks2 = [];
                for j = these_pks1(:)'
                    [~, vIdx] = ...
                        findDistToClosestFarBead(data, peaks.L6.beads{j}, peaks.L6.indices(j));
                    extra_pkbds = sort([peaks.L6.beads{j}(:); voxelToBeadMap.getBeadIdx(vIdx)]);
                    if fastIntersect(extra_pkbds, ridges1D.beads{i}, 'all')
                        these_pks2 = [these_pks2; j];
                    end
                end

                % remove peak that is already associated
                pk_inds = arrayfun(@(x) find(peaks.L6.indices == ridges1D.rp_key(i, x)), ...
                                   1 : ridges1D.rp_counter(i));
                these_pks2(fastIntersect(these_pks2, pk_inds, 'indices')) = [];

                % check if any peaks fit criteria; if not, remove ridge
                if isempty(these_pks2)
                    % remove ridge instead of creating new peak
                    remove_these_r1D(i) = true;
                % else, determine whether ridge contains enough voxels
                elseif length(ridges1D.indices{i}) <= 1
                    if length(these_pks2) > num_miss
                        % choose the peaks closest to the ridge voxel
                        pk_dists = arrayfun(@(x) norm(voxels(peaks.L6.indices(these_pks2(x)), :) ...
                                                      - ridges1D.indices{1}), ...
                                            1 : length(these_pks2));
                        [~, sort_ind] = sort(pk_dists);
                        for j = 1 : num_miss
                            pks2add(j) = peaks.L6.indices(these_pks2(sort_ind(j)));
                        end
                    elseif length(these_pks2) == num_miss
                        for j = 1 : num_miss
                            pks2add(j) = peaks.L6.indices(these_pks2(j));
                        end
                    else
                        % remove ridge instead of creating new peak
                        remove_these_r1D(i) = true;
                    end
                elseif num_miss == 1
                    if isempty(ridges1D.split{2, i}) || isempty(ridges1D.split{3, i})
                        % ridge does not have a bisecting normal plane;
                        % the following is a method for attempting to
                        % accurately find the other endpoint of the ridge
                        r1D_vox_end1 = [ridges1D.split{1, i}; ridges1D.split{4, i}];
                        % brute force method for finding voxel on opposite
                        % end of ridge
                        r1D_voxdist = arrayfun(@(x) norm(voxels(r1D_vox_end1, :) - ...
                                                         voxels(ridges1D.indices{i}(x), :)), ...
                                               1 : length(ridges1D.indices{i}));
                        [~, max_ind] = max(r1D_voxdist);
                        % use this voxel as a marker for the other end of ridge
                        r1D_vox_opp  = ridges1D.indices{i}(max_ind);

                        % List distances of ridge voxels to next-nearest bead
                        [dist, closestVoxelIdxs] = ...
                            findDistToClosestFarBead(data, ridges1D.beads{i}, ridges1D.indices{i});

                        % find the first voxel starting from smallest distance that lies on
                        % opposite end of ridge;
                        % remember row 'num_r1Dbds' in cumsum_mat actually
                        % refers to cumulative distance to bead 'num_r1Dbds + 1'
                        min_ind_opp = [];
                        % we're assuming that the voxel that is closest to
                        % its next-nearest bead is behaving the most like a
                        % peak
                        [~, sort_ind] = sort(dist); % sorting the cumsum distance of each voxel to its next-nearest bead
                        cntr = 1;
                        while isempty(min_ind_opp) && cntr <= length(sort_ind)
                            dist1 = norm(voxels(ridges1D.indices{i}(sort_ind(cntr)), :) - ...
                                         voxels(r1D_vox_end1, :));
                            dist2 = norm(voxels(ridges1D.indices{i}(sort_ind(cntr)), :) - ...
                                         voxels(r1D_vox_opp, :));
                            if dist2 < dist1
                                % voxel is closer to opposite end from current peak
                                min_ind_opp = sort_ind(cntr);
                            end
                            cntr = cntr + 1;
                        end
                    else
                        % use bisecting normal plane data
                        % List distances of ridge voxels to next-nearest bead
                        [dist, closestVoxelIdxs] = ...
                            findDistToClosestFarBead(data, ridges1D.beads{i}, ridges1D.indices{i});

                        % find the first voxel starting from smallest distance that lies on
                        % opposite end of ridge;
                        % remember row 'num_r1Dbds' in cumsum_mat actually
                        % refers to cumulative distance to bead 'num_r1Dbds + 1'
                        min_ind_opp = [];
                        [~, sort_ind] = sort(dist);
                        if isempty(ridges1D.split{1, i})
                            this_end = 2;
                        else
                            this_end = 3;
                        end
                        cntr = 1;
                        while isempty(min_ind_opp) && cntr <= length(sort_ind)
                            if binarySearch(ridges1D.split{this_end, i}, ...
                                            ridges1D.indices{i}(sort_ind(cntr))) > 0
                                % voxel lies on opposite end from current peak
                                min_ind_opp = sort_ind(cntr);
                            end
                            cntr = cntr + 1;
                        end
                    end

                    % use this voxel to match to peak
                    pot_beads = sort([ridges1D.beads{i}(:); ...
                                      voxelToBeadMap.getBeadIdx(closestVoxelIdxs(min_ind_opp))]);

                    % determine which potential peak has most overlapping beads
                    % with ridge
                    overlap = arrayfun(@(x) fastIntersect(peaks.L6.beads{these_pks2(x)}, ...
                                                          pot_beads, 'size'), ...
                                       1 : length(these_pks2));
                    pot_pks = these_pks2(overlap == max(overlap));
                    if isempty(pot_pks)
                        % remove ridge instead of creating new peak
                        remove_these_r1D(i) = true;
                    elseif length(pot_pks) == 1
                        pks2add = peaks.L6.indices(pot_pks);
                    else
                        % check whether potential peaks share a ridge;
                        % if not, combine peaks (this accommodates for
                        % dx issues with ridges2D). Otherwise choose
                        % peak that is closest to voxel
                        [og_pks, pot_pks, peaks.L6, peaks.pr1D_key, peaks.pr1D_counter, ...
                         ridges1D.indices, ridges1D.rp_key, ridges1D.split] = ...
                         combinePks(i, pot_pks, peaks.L6, ridges1D, ...
                                    peaks.pr1D_key, peaks.pr1D_counter, voxels);

                        % adjust peak indices in variable
                        pot_pks = fastIntersect(og_pks, pot_pks, 'indices');
                        % continue
                        if length(pot_pks) > 1
                            % choose peak that is closest to voxel
                            min_dist_ind = inf;
                            for j = pot_pks(:)'
                                dist = norm(voxels(peaks.L6.indices(j), :) - ...
                                            voxels(ridges1D.indices{i}(min_ind_opp), :));
                                if dist < min_dist_ind
                                    min_dist_ind = dist;
                                    pks2add = peaks.L6.indices(j);
                                end
                            end
                        else
                            pks2add = peaks.L6.indices(pot_pks);
                        end
                    end
                    % Check whether ridge already exists that flanks
                    % these peaks; if so, remove ridge
                    pks_pair = sort([pks2add, ridges1D.rp_key(i, 1)]);
                    if sum(arrayfun(@(x) fastIntersect(ridges1D.rp_key(x, 1 : 2), ...
                                         pks_pair, 'all'), ...
                                    ridges1D.pks2)) > 0
                        remove_these_r1D(i) = true;
                    end
                elseif num_miss == 2
                    min_ind = zeros(2, 1);
                    if isempty(ridges1D.split{2, i}) || isempty(ridges1D.split{3, i})
                        % ridge does not have a bisecting normal plane;
                        % List distances of ridge voxels to next-nearest bead
                        [dist, closestVoxelIdxs] = ...
                            findDistToClosestFarBead(data, ridges1D.beads{i}, ridges1D.indices{i});

                        % find the first voxel starting from smallest distance;
                        % remember row 'num_r1Dbds' in cumsum_mat actually
                        % refers to cumulative distance to bead 'num_r1Dbds + 1'
                        [~, sort_ind] = sort(dist);
                        min_ind(1) = sort_ind(1);
                        % brute force method for finding voxel on opposite
                        % end of ridge
                        r1D_voxdist = arrayfun(@(x) norm(voxels(ridges1D.indices{i}(min_ind(1)), :) - ...
                                                         voxels(ridges1D.indices{i}(x), :)), ...
                                               1 : length(ridges1D.indices{i}));
                        [~, max_ind] = max(r1D_voxdist);
                        % use this voxel as a marker for the other end of ridge
                        r1D_vox_opp  = ridges1D.indices{i}(max_ind);

                        cntr = 2;
                        while min_ind(2) == 0 && cntr <= length(sort_ind)
                            dist1 = norm(voxels(ridges1D.indices{i}(sort_ind(cntr)), :) - ...
                                         voxels(ridges1D.indices{i}(min_ind(1)), :));
                            dist2 = norm(voxels(ridges1D.indices{i}(sort_ind(cntr)), :) - ...
                                         voxels(r1D_vox_opp, :));
                            if dist2 < dist1
                                % voxel is closer to opposite end from current peak
                                min_ind(2) = sort_ind(cntr);
                            end
                            cntr = cntr + 1;
                        end
                    else
                        % use bisecting normal plane data
                        % List distances of ridge voxels to next-nearest bead
                        [dist, closestVoxelIdxs] = ...
                            findDistToClosestFarBead(data, ridges1D.beads{i}, ridges1D.indices{i});

                        % find the first voxel starting from smallest distance;
                        % remember row 'num_r1Dbds' in cumsum_mat actually
                        % refers to cumulative distance to bead 'num_r1Dbds + 1'
                        [~, sort_ind] = sort(dist);
                        % check which half of ridge the first voxel lies on
                        if binarySearch(ridges1D.split{2, i}, ...
                                        ridges1D.indices{i}(sort_ind(1))) > 0
                            min_ind(1) = sort_ind(1);
                            this_end = 3;
                            ii = 2;
                        else
                            min_ind(2) = sort_ind(1);
                            this_end = 2;
                            ii = 1;
                        end

                        cntr = 2;
                        while min_ind(ii) == 0 && cntr <= length(sort_ind)
                            if binarySearch(ridges1D.split{this_end, i}, ...
                                            ridges1D.indices{i}(sort_ind(cntr))) > 0
                                % voxel lies on opposite end from current peak
                                min_ind(ii) = sort_ind(cntr);
                            end
                            cntr = cntr + 1;
                        end
                    end
                    % use these voxels to match to peak
                    for j = 1 : 2
                        pot_beads = sort([ridges1D.beads{i}(:); ...
                                          voxelToBeadMap.getBeadIdx(closestVoxelIdxs(min_ind(j)))]);

                        % determine which potential peak has most overlapping beads
                        % with ridge
                        overlap = arrayfun(@(x) fastIntersect(peaks.L6.beads{these_pks2(x)}, ...
                                                              pot_beads, 'size'), ...
                                           1 : length(these_pks2));
                        pot_pks = these_pks2(overlap == max(overlap));
                        if isempty(pot_pks)
                            % remove ridge instead of creating new peak
                            remove_these_r1D(i) = true;
                        elseif length(pot_pks) == 1
                            pks2add(j) = peaks.L6.indices(pot_pks);
                        else
                            % check whether potential peaks share a ridge;
                            % if not, combine peaks (this accommodates for
                            % dx issues with ridges2D). Otherwise choose
                            % peak that is closest to voxel
                            [og_pks, pot_pks, peaks.L6, peaks.pr1D_key, peaks.pr1D_counter, ...
                             ridges1D.indices, ridges1D.rp_key, ridges1D.split] = ...
                             combinePks(i, pot_pks, peaks.L6, ridges1D, ...
                                        peaks.pr1D_key, peaks.pr1D_counter, voxels);

                            % adjust peak indices in variable
                            these_pks2 = fastIntersect(og_pks, these_pks2, 'indices');
                            pot_pks    = fastIntersect(og_pks, pot_pks, 'indices');
                            % continue
                            if length(pot_pks) > 1
                                % choose peak that is closest to voxel
                                min_dist_ind = inf;
                                for k = pot_pks(:)'
                                    dist = norm(voxels(peaks.L6.indices(k), :) - ...
                                                voxels(ridges1D.indices{i}(min_ind(j)), :));
                                    if dist < min_dist_ind
                                        min_dist_ind = dist;
                                        pks2add(j) = peaks.L6.indices(k);
                                    end
                                end
                            else
                                pks2add(j) = peaks.L6.indices(pot_pks);
                            end
                        end
                    end
                    % Check whether ridge already exists that flanks
                    % these peaks; if so, remove ridge
                    pks_pair = sort(pks2add);
                    if sum(arrayfun(@(x) fastIntersect(ridges1D.rp_key(x, 1 : 2), ...
                                         pks_pair, 'all'), ...
                                    ridges1D.pks2)) > 0
                        remove_these_r1D(i) = true;
                    end
                elseif num_miss > 2
                    error('There are incorrectly more than 2 peaks missing for a ridge.')
                end
                % Update ridge
                if length(pks2add) == 1 && pks2add ~= 0
                    % split
                    if isempty(ridges1D.split{1, i})
                        ridges1D.split{1, i} = pks2add;
                    else
                        ridges1D.split{4, i} = pks2add;
                    end
                    % key & counter
                    ridges1D.rp_key(i, 1 : 2) = sort([ridges1D.rp_key(i, 1), ...
                                                      pks2add]);
                    ridges1D.rp_counter(i) = 2;
                elseif length(pks2add) == 2 && sum(pks2add) ~= 0
                    % split
                    ridges1D.split{1, i} = pks2add(1);
                    ridges1D.split{4, i} = pks2add(2);
                    % key & counter
                    ridges1D.rp_key(i, 1 : 2) = sort(pks2add)';
                    ridges1D.rp_counter(i) = 2;
                end
                % Update ridges1D types
                ridges1D.pks0 = find(ridges1D.rp_counter == 0);
                ridges1D.pks1 = find(ridges1D.rp_counter == 1);
                ridges1D.pks2 = find(ridges1D.rp_counter >= 2);
            end
        end
        % Remove necessary ridges
        ridges1D.num = ridges1D.num - sum(remove_these_r1D);
        ridges1D.indices(remove_these_r1D)       = [];
        ridges1D.beads(remove_these_r1D)         = [];
        ridges1D.rp_key(remove_these_r1D, :)     = [];
        ridges1D.rp_counter(remove_these_r1D)    = [];
        ridges1D.doors(remove_these_r1D)         = [];
        ridges1D.lowest_avg(remove_these_r1D, :) = [];
        ridges1D.min(remove_these_r1D)           = [];
        ridges1D.bead_voxs(remove_these_r1D, :)  = [];
        ridges1D.split(:, remove_these_r1D)      = [];
        ridges1D.edge_bool(remove_these_r1D)     = [];

        % recreate peaks.pr1D_key
        peaks.pr1D_key       = zeros(length(peaks.L6.indices), 16);
        peaks.pr1D_key(:, 1) = peaks.L6.indices;
        peaks.pr1D_counter   = ones(size(peaks.pr1D_key, 1), 1);
        % convert ridges1D.rp_key to peak indices
        r1D_pkvoxs         = ridges1D.rp_key(:, 1 : 2);
        ridges1D.rp_keyind = cell2mat(arrayfun(@(x) ...
                                      [binarySearch(peaks.L6.indices, r1D_pkvoxs(x, 1)), ...
                                       binarySearch(peaks.L6.indices, r1D_pkvoxs(x, 2))], ...
                                      1 : size(r1D_pkvoxs, 1), 'UniformOutput', false)');
        % fill key
        for i = 1 : ridges1D.num
            for j = 1 : min(ridges1D.rp_counter(i), 2)
                if ridges1D.rp_keyind(i, j) ~= -1
                    this_pkrow = ridges1D.rp_keyind(i, j);
                    peaks.pr1D_counter(this_pkrow) = peaks.pr1D_counter(this_pkrow) + 1;
                    peaks.pr1D_key(this_pkrow, peaks.pr1D_counter(this_pkrow)) = i; % already sorted
                end
            end
        end
        peaks.pr1D_counter = peaks.pr1D_counter - 1;

        % Update ridges1D types
        ridges1D.pks0 = find(ridges1D.rp_counter == 0);
        ridges1D.pks1 = find(ridges1D.rp_counter == 1);
        ridges1D.pks2 = find(ridges1D.rp_counter >= 2);
    end


    % convert ridges1D.rp_key to peak indices
    r1D_pkvoxs         = ridges1D.rp_key(:, 1 : 2);
    ridges1D.rp_keyind = cell2mat(arrayfun(@(x) ...
                                  [binarySearch(peaks.L6.indices, r1D_pkvoxs(x, 1)), ...
                                   binarySearch(peaks.L6.indices, r1D_pkvoxs(x, 2))], ...
                                  1 : size(r1D_pkvoxs, 1), 'UniformOutput', false)');
    % store edge ridges
    ridges1D.edge = find(ridges1D.edge_bool);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Force peaks:');

    time_log(timeLogIdx).Name = 'Force peaks';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % COMBINE PEAKS
    % For peaks that are surrounding a non-existent 'true' peak (i.e., peak
    % doesn't lie on the grid), choose the highest existing local peak
    tStart = tic;

    % Next peaks level:
    peaks.L7 = peaks.L6;

    % we attempt to identify peaks surrounding true peaks by looking for
    % peaks that are associated with < 3 ridges1D
    peaks_2r1D = find(peaks.pr1D_counter < 3);
    % cluster peaks that share neighboring beads
    queue = peaks_2r1D(:)';
    clustered_pks = cell(length(queue), 1);
    cntr = 1;
    while ~isempty(queue)
        lead = queue(1);
        clustered_pks{cntr} = lead;
        queue(1) = [];
        remove_pos = false(length(queue), 1);
        max_size = peaks.L6.numBeads(lead);
        for i = 1 : length(queue)
            % check if share 3 out of 4 or 4 out of 5, etc. neighboring beads
            if fastIntersect(peaks.L6.beads{lead}, peaks.L6.beads{queue(i)}, 'size') >= max_size - 1
                % if so, ensure peaks don't already share ridge1D
                if ~fastIntersect(peaks.pr1D_key(lead, 2 : peaks.pr1D_counter(lead) + 1), ...
                                  peaks.pr1D_key(queue(i), 2 : peaks.pr1D_counter(queue(i)) + 1), 'bool')
                    clustered_pks{cntr} = [clustered_pks{cntr}; queue(i)];
                    remove_pos(i) = true;
                end
            end
        end
        queue(remove_pos) = [];
        cntr = cntr + 1;
    end
    clustered_pks = clustered_pks(~cellfun('isempty', clustered_pks));
    % do not include peaks that are not part of a cluster
    clustered_pks = clustered_pks(cellfun(@(x) length(x) > 1, clustered_pks));
    for i = 1 : length(clustered_pks)
        % select highest peak
        [~, ind]     = max(data.EDT(peaks.L6.indices(clustered_pks{i})));
        chosen_pkrow = clustered_pks{i}(ind);
        % update data to remove other peaks
        peaks.L7.beads{chosen_pkrow}    = unique(cell2mat(peaks.L7.beads(clustered_pks{i})));
        peaks.L7.numBeads(chosen_pkrow) = length(peaks.L7.beads{chosen_pkrow});
        % collect total number of ridges1D and update key
        updated_ridges = nonzeros(unique(peaks.pr1D_key(clustered_pks{i}, 2 : end)))';
        peaks.pr1D_key(chosen_pkrow, 2 : length(updated_ridges) + 1) = updated_ridges;
        peaks.pr1D_counter(chosen_pkrow) = length(updated_ridges);
        % store only old peaks
        clustered_pks{i}(ind) = [];
        % update ridge1D data
        for j = updated_ridges
            og_ridges    = ridges1D.rp_key(j, 1 : ridges1D.rp_counter(j));
            for k = clustered_pks{i}(:)'
                og_ridges(og_ridges == peaks.L7.indices(k)) = peaks.L7.indices(chosen_pkrow);
            end
            og_ridges = unique(og_ridges(:))';
            ridges1D.rp_counter(j)   = length(og_ridges);
            ridges1D.rp_key(j, :)    = 0 * ridges1D.rp_key(j, :); % reset
            ridges1D.rp_key(j, 1 : ridges1D.rp_counter(j)) = og_ridges;
            for k = [1, 4]
                if binarySearch(peaks.L7.indices(clustered_pks{i}), ridges1D.split{k, j}) > 0
                    ridges1D.split{k, j} = peaks.L7.indices(chosen_pkrow);
                end
            end
        end
    end
    % update peak number
    peaks.L7.num = peaks.L7.num - length(cell2mat(clustered_pks));
    % remove old peaks
    remove_these_pks = cell2mat(clustered_pks);
    peaks.L7.indices(remove_these_pks)  = [];
    peaks.L7.beads(remove_these_pks)     = [];
    peaks.L7.numBeads(remove_these_pks)  = [];
    peaks.L7.EDT(remove_these_pks)       = [];
    peaks.pr1D_key(remove_these_pks, :)  = [];
    peaks.pr1D_counter(remove_these_pks) = [];

    % remove ridges associated with duplicate peaks
    rmv_duplicate = [];
    for i = 1 : size(ridges1D.rp_key, 1)
        % not even considering 3 duplicates...
        if (ridges1D.rp_key(i, 1) == ridges1D.rp_key(i, 2)) && (ridges1D.rp_key(i, 1) ~= 0)
            rmv_duplicate = [rmv_duplicate; i];
        end
    end
    ridges1D.rp_key(rmv_duplicate, 2)    = 0;
    ridges1D.rp_counter(rmv_duplicate)   = 1;
    ridges1D.rp_keyind(rmv_duplicate, 2) = 0;
    ridges1D.split{4, i} = []; % technically should check which end the peak is closest to

    % Update ridge1D types
    ridges1D.pks0 = find(ridges1D.rp_counter == 0);
    ridges1D.pks1 = find(ridges1D.rp_counter == 1);
    ridges1D.pks2 = find(ridges1D.rp_counter >= 2);

    % convert ridges1D.rp_key to peak indices
    r1D_pkvoxs         = ridges1D.rp_key(:, 1 : 2);
    ridges1D.rp_keyind = cell2mat(arrayfun(@(x) ...
                                  [binarySearch(peaks.L7.indices, r1D_pkvoxs(x, 1)), ...
                                   binarySearch(peaks.L7.indices, r1D_pkvoxs(x, 2))], ...
                                  1 : size(r1D_pkvoxs, 1), 'UniformOutput', false)');

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Combine peaks:');

    time_log(timeLogIdx).Name = 'Combine peaks';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    %%%%%%%************************************************%%%%%%%
    %%%%%***************** CONNECTING PEAKS *****************%%%%%
    %%%%%%%************************************************%%%%%%%

    % Determine which peaks to combine based on two methods:
    % Method #1: check whether EDT of min point along ridge is > 80% of avg
    %            EDT between both peaks
    % Method #2: check if physical distance between peaks is < sum of each
    %            EDT (i.e., radii)

    tStart = tic;

    rp = ridges1D.rp_key(ridges1D.pks2, 1 : 2);
    pks_graph = ridges1D.rp_keyind(ridges1D.pks2, :);

%     % Method #1
%     EDT_cutoff = dip_percent * mean(peaks.L7.EDT(pks_graph), 2);
%     r1D_min    = edt_full(ridges1D.min(ridges1D.pks2));
%     ridges1D.threshBool1 = r1D_min > EDT_cutoff; % bool only applies to ridges1D.pks2 vector!

    % Method #2
    phys_dist = vecnorm(voxels(rp(:, 1), :) - voxels(rp(:, 2), :), 2, 2);
    EDT_radii = arrayfun(@(x) peaks.L7.EDT(pks_graph(x, 1)) + ...
                              peaks.L7.EDT(pks_graph(x, 2)), ...
                              1 : size(pks_graph, 1))';
    ridges1D.threshBool2 = phys_dist < EDT_radii; % bool only applies to ridges1D.pks2 vector!
    % For forcing no connections:
    if dip_percent >= 1
        ridges1D.threshBool2 = false(length(ridges1D.threshBool2), 1);
    end

    % store ridges that connect peaks
    ridges1D.connected = ridges1D.pks2(ridges1D.threshBool2);

    % create graph for connected peaks
    pks_graph_conn = pks_graph(ridges1D.threshBool2, :); % has length of ridge1D.pk2

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Threshold to combine peaks:');

    time_log(timeLogIdx).Name = 'Threshold to combine peaks';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    tStart = tic;
    % Determine peaks per subunit using adjacency matrix for connected peaks
    n = max(pks_graph_conn(:));
    peaks.adjacency = sparse(pks_graph_conn(:, 1), pks_graph_conn(:, 2), ...
                             ones(size(pks_graph_conn, 1), 1), ...
                             n, n);

    % Store data in struct
    subs_clust           = struct;
    subs_clust.peaks.row = cell(peaks.L7.num, 1);
    subs_clust.peaks.ind = cell(peaks.L7.num, 1);
    visit_bool           = false(peaks.L7.num, 1);
    clust_c              = 1;

    while ~all(visit_bool)
        current_cache                 = find(~visit_bool, true, 'first');
        subs_clust.peaks.row{clust_c} = current_cache(1);
        visit_bool(current_cache(1))  = true;
        while ~isempty(current_cache)
            peak_current     = current_cache(1);
            current_cache(1) = [];
            % find other peaks connecting to current peak (sharing ridges1D)
            neigh_peaks_row = find(peaks.adjacency(peak_current, :));
            neigh_peaks_col = find(peaks.adjacency(:, peak_current));
            neigh_peaks = [neigh_peaks_row(:); neigh_peaks_col(:)];
            % remove peaks that have already been visited
            neigh_peaks = neigh_peaks(~visit_bool(neigh_peaks));
            visit_bool(neigh_peaks) = true;
            % add to current cluster
            subs_clust.peaks.row{clust_c} = [subs_clust.peaks.row{clust_c}; neigh_peaks(:)];
            current_cache = [current_cache; neigh_peaks(:)];
        end
        subs_clust.peaks.row{clust_c} = sort(subs_clust.peaks.row{clust_c});
        subs_clust.peaks.ind{clust_c} = peaks.L7.indices(subs_clust.peaks.row{clust_c});
        clust_c = clust_c + 1;
    end
    % remove empty rows
    subs_clust.peaks.row = subs_clust.peaks.row(~cellfun('isempty', subs_clust.peaks.row));
    subs_clust.peaks.ind = subs_clust.peaks.ind(1 : length(subs_clust.peaks.row));
    % store number of subunits
    num_subs = length(subs_clust.peaks.row);

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Cluster peaks per subunit:');

    time_log(timeLogIdx).Name = 'Cluster peaks per subunit';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    %%%%%%%***********************************************%%%%%%%
    %%%%%***************** RIDGE1D LENGTHS *****************%%%%%
    %%%%%%%***********************************************%%%%%%%

    tStart = tic;

    ridges1D.lengths = zeros(ridges1D.num, 1);
    % Compute length for each 1D ridge
    for i = 1 : ridges1D.num
        if ridges1D.rp_counter(i) == 0
            % no peak flanks ridge, so choose an 'end point'
            r1D_edgeind = fastIntersect(ridges1D.indices{i}, edge_ind, 'elements');
            if ~isempty(r1D_edgeind)
                % choose point that intersects with edge of void space
                this_r1Dind = r1D_edgeind(1);
            else
                % if no indices intersect with edge of void space, just choose
                % first index in list of points
                this_r1Dind = ridges1D.indices{i}(1);
            end
            ridges1D.lengths(i) = ridgeLength(ridges1D.indices{i}, this_r1Dind, voxels, dx);
        else
            try
                ridges1D.lengths(i) = ridgeLength(ridges1D.indices{i}, ridges1D.rp_key(i, 1 : ridges1D.rp_counter(i)), voxels, dx);
            catch
                continue;
            end
        end
    end
    % Store peaks that share 1D ridges and include ridge lengths (for
    % creating graph object)
    ridge_length_graph   = ridges1D.lengths(ridges1D.pks2);
    peaks.pks_graph_length = [pks_graph, ridge_length_graph];

    % Get path data for full domain
    [center_pk, path_nodes, path_length, path_r1Ds, path_tortuosity, path_necks, path_doors, path_r1D_doors] = ...
                pathsCenter(peaks.L7.indices, peaks.pks_graph_length, ...
                                ridges1D, edge_ind, voxels, domain);

    data.paths.num_paths    = numel(path_nodes);
    data.paths.center_pk    = center_pk;
    data.paths.node_pks     = path_nodes;
    data.paths.path_lengths = path_length;
    data.paths.path_r1Ds    = path_r1Ds;
    data.paths.tortuosity   = path_tortuosity;
    data.paths.neck_data    = path_necks;
    data.paths.door_data    = path_doors;
    data.paths.r1D_w_doors  = path_r1D_doors;

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Lengths of ridges1D:');

    time_log(timeLogIdx).Name = 'Lengths of ridges1D';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    %%%%%%%********************************************************%%%%%%%
    %%%%%***************** RIDGE1D DATA PER SUBUNIT *****************%%%%%
    %%%%%%%********************************************************%%%%%%%

    % Determine ridges1D per subunit using peaks data
    subs_clust.ridges1D = cell(num_subs, 1);
    subs_clust.ridges1D_lengths = cell(num_subs, 1);
    for i = 1 : num_subs
        subs_clust.ridges1D{i} = unique(cell2mat(arrayfun(@(x) ...
                                                 peaks.pr1D_key(x, 2 : 1 + peaks.pr1D_counter(x))', ...
                                                 subs_clust.peaks.row{i}, 'UniformOutput', false)));
        % store ridge1D # and associated length
        subs_clust.ridges1D_lengths{i} = zeros(length(subs_clust.ridges1D{i}), 2);
        subs_clust.ridges1D_lengths{i}(:, 1) = subs_clust.ridges1D{i};
    end

    tStart = tic;
    % Determine skeleton (or backbone) voxels per subunit and track neighbors
    subs_clust.skeleton    = cell(num_subs, 1);
    % Gather 'center' of subunit, which entails peaks and connecting (full) ridges1D
    subs_clust.center      = cell(num_subs, 1);
    subs_clust.center_r1Dmin_diams = cell(num_subs, 1);
    % Store ridges containing doors
    subs_clust.door_ridges = cell(num_subs, 1);
    subs_clust.edge_ridges = cell(num_subs, 1);
    subs_clust.hall_ridges = cell(num_subs, 1);
    for i = 1 : num_subs
        fullhalf_bool = fastIntersect(subs_clust.ridges1D{i}, ridges1D.connected, 'boolvec');
        edge_bool = fastIntersect(subs_clust.ridges1D{i}, ridges1D.edge, 'boolvec');
        % these are the 'half ridges' associated with the subunit
        get_half = subs_clust.ridges1D{i}(~fullhalf_bool & ~edge_bool);

        % fill neighbors graph data
        subs_clust.hall_ridges{i} = get_half(:);
        % Add voxels for HALF ridges
        if ~isempty(get_half)
            for j = get_half(:)'
                % check top half for corresponding peak within subunit
                if fastIntersect(ridges1D.split{1, j}, subs_clust.peaks.ind{i}, 'bool')
                    subs_clust.skeleton{i} = [subs_clust.skeleton{i}; ridges1D.split{2, j}];
                    % compute length and store in matrix
                    this_ridge = binarySearch(subs_clust.ridges1D_lengths{i}(:, 1), j);
                    % add ridge1D min point as other anchor point
                    if ~isempty(ridges1D.split{2, j})
                        subs_clust.ridges1D_lengths{i}(this_ridge, 2) = ...
                                ridgeLength(ridges1D.split{2, j}, [ridges1D.split{1, j}, ridges1D.min(j)], voxels, dx);
                    else
                        subs_clust.ridges1D_lengths{i}(this_ridge, 2) = 0;
                    end
                end
                % check bottom half for corresponding peak within subunit
                if fastIntersect(ridges1D.split{4, j}, subs_clust.peaks.ind{i}, 'bool')
                    subs_clust.skeleton{i} = [subs_clust.skeleton{i}; ridges1D.split{3, j}];
                    % compute length and store in matrix
                    this_ridge = binarySearch(subs_clust.ridges1D_lengths{i}(:, 1), j);
                    % add ridge1D min point as other anchor point
                    if ~isempty(ridges1D.split{3, j})
                        subs_clust.ridges1D_lengths{i}(this_ridge, 2) = ...
                                ridgeLength(ridges1D.split{3, j}, [ridges1D.split{4, j}, ridges1D.min(j)], voxels, dx);
                    else
                        subs_clust.ridges1D_lengths{i}(this_ridge, 2) = 0;
                    end
                end
            end
        end
        % Add voxels for EDGE ridges
        get_edge = subs_clust.ridges1D{i}(edge_bool);
        edge_voxs = cell2mat(ridges1D.indices(get_edge));
        subs_clust.edge_ridges{i} = get_edge;
        subs_clust.skeleton{i} = [subs_clust.skeleton{i}; edge_voxs(:)];
        for j = get_edge(:)'
            % compute length and store in matrix
            this_ridge = binarySearch(subs_clust.ridges1D_lengths{i}(:, 1), j);
            if ~isempty(ridges1D.indices{j})
                subs_clust.ridges1D_lengths{i}(this_ridge, 2) = ... % edge ridges should contain single peak in 1st column of rp_key
                        ridgeLength(ridges1D.indices{j}, ridges1D.rp_key(j, 1), voxels, dx);
            else
                subs_clust.ridges1D_lengths{i}(this_ridge, 2) = 0;
            end
        end
        % Add voxels for FULL ridges
        get_full = subs_clust.ridges1D{i}(fullhalf_bool);
        full_voxs = cell2mat(ridges1D.indices(get_full));
        subs_clust.skeleton{i} = [subs_clust.skeleton{i}; full_voxs(:)];
        subs_clust.center{i}   = [subs_clust.center{i}; full_voxs(:)];
        for j = get_full(:)'
            % compute length and store in matrix
            this_ridge = binarySearch(subs_clust.ridges1D_lengths{i}(:, 1), j);
            if ~isempty(ridges1D.indices{j})
                subs_clust.ridges1D_lengths{i}(this_ridge, 2) = ... % edge ridges should contain single peak in 1st column of rp_key
                        ridgeLength(ridges1D.indices{j}, ridges1D.rp_key(j, 1:2)', voxels, dx);
            else
                subs_clust.ridges1D_lengths{i}(this_ridge, 2) = norm(voxels(peaks.L7.indices(ridges1D.rp_key(j, 1)), :) - ...
                                                                     voxels(peaks.L7.indices(ridges1D.rp_key(j, 2)), :));
            end
            % store diameter (EDT*2) at minimum (i.e., at door) along full ridge
            subs_clust.center_r1Dmin_diams{i} = [subs_clust.center_r1Dmin_diams{i}; ...
                                                 ridges1D.doors{j}.radius * 2];
        end
        % Add voxels for PEAKS
        subs_clust.skeleton{i} = [subs_clust.skeleton{i}; subs_clust.peaks.ind{i}];
        subs_clust.center{i}   = [subs_clust.center{i}; subs_clust.peaks.ind{i}];
        % sort
        subs_clust.skeleton{i} = sort(subs_clust.skeleton{i});
        subs_clust.center{i}   = sort(subs_clust.center{i});
        % Gather ridges1D that contain interior doors of subunit
        subs_clust.door_ridges{i} = get_half;
    end

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Gather subunit skeleton:');

    time_log(timeLogIdx).Name = 'Gather subunit skeleton';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    % Subunit bead adjacency matrix for associating ridges2D to subs
    neighbors_beads     = cell(num_subs, 1);
    subs_adjacencybeads = cell(num_subs, 1);
    for i = 1 : num_subs
        neighbors_beads{i} = unique(cell2mat(arrayfun(@(x) peaks.L7.beads{x}, ...
                                                                    subs_clust.peaks.row{i}, ...
                                                                    'UniformOutput', false)));
        % creating 2-column array of: subunit number (col 1) & bead number (col 2)
        % for later adjacency matrix
        subs_adjacencybeads{i} = [i * ones(length(neighbors_beads{i}), 1), ...
                                  neighbors_beads{i}];
    end

    tStart = tic;
    % Sparse subunit-beads adjacency matrix
    subs_adjacencybeads = cell2mat(subs_adjacencybeads);
    subs_adjacencybeads = sparse(subs_adjacencybeads(:, 1), subs_adjacencybeads(:, 2), ...
                                 ones(size(subs_adjacencybeads, 1), 1), ...
                                 num_subs, num_beads);

    % Determine ridges2D per subunit using subunit-beads adjacency matrix
    subs_clust.ridges2D = cell(num_subs, 1);
    for i = 1 : size(ridges2D.beads, 1)
        subs1 = find(subs_adjacencybeads(:, ridges2D.beads{i}(1)));
        subs2 = find(subs_adjacencybeads(:, ridges2D.beads{i}(2)));
        % intersect to find all subunits that contain ridge2D #i
        these_subs = fastIntersect(subs1, subs2, 'elements');
        filling_key = arrayfun(@(x) vertcat(subs_clust.ridges2D{x}, i), ...
                                            these_subs, 'UniformOutput', false);
        subs_clust.ridges2D(these_subs) = filling_key;
    end

    tElapsed = toc(tStart);
    % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Associate ridges2D to subs:');

    time_log(timeLogIdx).Name = 'Associate ridges2D to subs';
    time_log(timeLogIdx).Time = tElapsed;
    time_log(timeLogIdx).Units = 'sec';
    printTimeInfo(time_log(timeLogIdx));
    timeLogIdx = timeLogIdx + 1;

    %%%%%%%************************************************%%%%%%%
    %%%%%***************** FORMING SUBUNITS *****************%%%%%
    %%%%%%%************************************************%%%%%%%

    if ~isscalarzero(peaks.L7.indices) && ~isempty(num_subs)
        tStart = tic;
        % remove skeleton voxels from void space
        vsVoxels_noSkel = vsVoxels;
        skeleton_cum = cell2mat(arrayfun(@(x) [subs_clust.skeleton{x}, ...
                                               x * ones(length(subs_clust.skeleton{x}), 1)]', ...
                                               1 : length(subs_clust.skeleton), ...
                                               'UniformOutput', false))';
        noskel_bool = fastIntersect(vsVoxels_noSkel, sort(skeleton_cum(:, 1)), 'boolvec');
        vsVoxels_noSkel = vsVoxels_noSkel(~noskel_bool);

        % For every void space voxel, determine the nearest skeleton
        knnsearch_ptCloud = knnsearch(voxels(skeleton_cum(:, 1), :), ...
                                      voxels(vsVoxels_noSkel, :), ...
                                      'k', 1, 'NSMethod', 'kdtree');

        % convert knnsearch output to corresponding subunit number
        knnsearch_convert = skeleton_cum(knnsearch_ptCloud, 2);
        % collect complete subunit voxels
        subs_clust.indices = subs_clust.skeleton;
        for i = 1 : num_subs
            subs_clust.indices{i} = sort([subs_clust.indices{i}; ...
                                          vsVoxels_noSkel(knnsearch_convert ==i)]);
        end

        % Get indices of 2D ridges within each subunit (used to compute local
        % thickness)
        subs_clust.skeleton2D = cell(num_subs, 1);
        for i = 1 : num_subs
            % Important to note that ridges2D.indices{i} is already sorted!
            tmp = arrayfun(@(i) ridges2D.indices{i}, subs_clust.ridges2D{i}, 'UniformOutput', false);
            tmp = cellfun(@(x) fastIntersect(x, subs_clust.indices{i}, 'elements'), tmp, 'UniformOutput', false);
            subs_clust.skeleton2D{i} = sort(vertcat(tmp{:}));
        end

        tElapsed = toc(tStart);
        % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Form subunits:');

        time_log(timeLogIdx).Name = 'Form subunits';
        time_log(timeLogIdx).Time = tElapsed;
        time_log(timeLogIdx).Units = 'sec';
        printTimeInfo(time_log(timeLogIdx));
        timeLogIdx = timeLogIdx + 1;

        % Get edge voxels of each subunit
        tStart = tic;
        subs_clust.edgeind = cell(num_subs, 1);
        subs_clust.edgeEDT = cell(num_subs, 1);
        edge_subs_log = false(num_subs, 1);
        edge_clust = cell(length(edgesubs2comb_2Dind), 1);
        no_comb_edgesubs = [];
        for i = 1 : num_subs
            [subs_clust.edgeind{i}, subs_clust.edgeEDT{i}] = getEdgeData(subs_clust.indices{i}, ...
                                                                         voxels, shape, dx, ...
                                                                         subs_clust.center{i});
            no_comb = 0;
            % edge subunit or not
            if fastIntersect(subs_clust.edgeind{i}, edge_ind, 'bool')
                edge_subs_log(i) = true;
                % check to which edge cluster sub belongs
                if combine_edge_subs
                    for j = 1 : length(edgesubs2comb_2Dind)
                        if fastIntersect(edgesubs2comb_2Dind{j}, subs_clust.edgeind{i}, 'bool')
                            edge_clust{j} = [edge_clust{j}; i];
                            no_comb = 1;
                        end
                    end
                end
                if no_comb == 0
                    no_comb_edgesubs = [no_comb_edgesubs; i];
                end
            else
                edge_subs_log(i) = false;
            end
        end
        % combine clusters that share subs
        edge_clust_comb = cell(length(edge_clust), 1);
        check_it        = true(length(edge_clust), 1);
        cntr = 1;
        for i = 1 : length(edge_clust)
            if check_it(i)
                current_group = edge_clust{i};
                trkr = i + 1;
                while trkr <= length(edge_clust)
                    if check_it(trkr) && fastIntersect(current_group, ...
                                                       edge_clust{trkr}, 'bool')
                        current_group = sort([current_group; edge_clust{trkr}]);
                        check_it(trkr) = false;
                        trkr = i + 1;
                    else
                        trkr = trkr + 1;
                    end
                end
                edge_clust_comb{cntr} = unique(current_group);
                cntr = cntr + 1;
            end
        end
        edge_clust_comb = edge_clust_comb(~cellfun('isempty', edge_clust_comb));
        edge_clust_all = [edge_clust_comb; mat2cell(no_comb_edgesubs, ones(numel(no_comb_edgesubs), 1))];

        % Match up edge peaks to edge subunits
        % get edge indices of edge subs
        edge_clust_inds = cell(numel(edge_clust_all), 1);
        edge_clust_pks  = cell(numel(edge_clust_all), 1);
        remove_edge_sub = false(numel(edge_clust_all), 1);
        for i = 1 : numel(edge_clust_all)
            % combine subunit indices if necessary
            if length(edge_clust_all{i}) > 1
                edge_clust_inds{i} = sort(cell2mat(subs_clust.indices(edge_clust_all{i})));
            else
                edge_clust_inds{i} = subs_clust.indices{edge_clust_all{i}};
            end
            % scan through peaks to find matches
            edge_peaks_cache = peaks_e.indices;
            remove_pks = false(numel(edge_peaks_cache), 1);
            % check each peak
            for j = 1 : numel(edge_peaks_cache)
                peak_found = binarySearch(edge_clust_inds{i}, edge_peaks_cache(j));
                if peak_found > 0
                    edge_clust_pks{i} = [edge_clust_pks{i}; edge_peaks_cache(j)];
                    remove_pks(j) = true;
                end
            end
            edge_peaks_cache(remove_pks) = [];
            edge_clust_pks{i} = sort(edge_clust_pks{i});
            if isempty(edge_clust_pks{i})
                % if there is no edge peak associated with edge subunit,
                % then remove subunit from edge subunit list
                remove_edge_sub(i) = true;
            end
        end
        edge_clust_all(remove_edge_sub) = [];
        edge_clust_inds(remove_edge_sub) = [];

        numEdgeSubs = length(edge_clust_all);

        % Cluster edge subunits that share the same open space
        subs_clust_e.peaks.row        = cell(numEdgeSubs, 1);
        subs_clust_e.peaks.ind        = cell(numEdgeSubs, 1);
        subs_clust_e.ridges1D         = cell(numEdgeSubs, 1);
        subs_clust_e.ridges1D_lengths = cell(numEdgeSubs, 1);
        subs_clust_e.skeleton         = cell(numEdgeSubs, 1);
        subs_clust_e.skeleton2D       = cell(numEdgeSubs, 1);
        subs_clust_e.center           = cell(numEdgeSubs, 1);
        subs_clust_e.center_r1Dmin_diams = cell(numEdgeSubs, 1);
        subs_clust_e.door_ridges      = cell(numEdgeSubs, 1);
        subs_clust_e.edge_ridges      = cell(numEdgeSubs, 1);
        subs_clust_e.hall_ridges      = cell(numEdgeSubs, 1);
        subs_clust_e.ridges2D         = cell(numEdgeSubs, 1);
        subs_clust_e.indices          = cell(numEdgeSubs, 1);
        subs_clust_e.edgeind          = cell(numEdgeSubs, 1);
        subs_clust_e.edgeEDT          = cell(numEdgeSubs, 1);
        for i = 1 : numEdgeSubs
            if length(edge_clust_all{i}) > 1
                % Combine everything
                subs_clust_e.peaks.row{i}    = unique(cell2mat(subs_clust.peaks.row(edge_clust_all{i})));
                subs_clust_e.peaks.ind{i}    = unique(cell2mat(subs_clust.peaks.ind(edge_clust_all{i})));
                subs_clust_e.ridges1D{i}     = unique(cell2mat(subs_clust.ridges1D(edge_clust_all{i})));
                % >> need to combine 2 columns for ridges1D_length data
                ae = sortrows(cell2mat(subs_clust.ridges1D_lengths(edge_clust_all{i})));
                subs_clust_e.ridges1D_lengths{i} = zeros(size(ae, 1), 2);
                tkr = 1;
                ptr = 1;
                while tkr <= size(ae, 1)
                    tkr_0 = tkr;
                    if tkr < size(ae, 1)
                        while ae(tkr, 1) == ae(tkr + 1, 1)
                            tkr = tkr + 1;
                            if tkr >= size(ae, 1)
                                tkr = tkr  -  1;
                                break;
                            end
                        end
                    end
                    subs_clust_e.ridges1D_lengths{i}(ptr, 1) = ae(tkr, 1);
                    subs_clust_e.ridges1D_lengths{i}(ptr, 2) = sum(ae(tkr_0 : tkr, 2));
                    tkr = tkr + 1;
                    ptr = ptr + 1;
                end
                subs_clust_e.ridges1D_lengths{i}(ptr : end, :) = [];
                subs_clust_e.skeleton{i}     = unique(cell2mat(subs_clust.skeleton(edge_clust_all{i})));
                subs_clust_e.skeleton2D{i}   = unique(cell2mat(subs_clust.skeleton2D(edge_clust_all{i})));
                subs_clust_e.center{i}       = unique(cell2mat(subs_clust.center(edge_clust_all{i})));
                subs_clust_e.center_r1Dmin_diams{i} = unique(cell2mat(subs_clust.center_r1Dmin_diams(edge_clust_all{i})));
                subs_clust_e.door_ridges{i}  = unique(cell2mat(subs_clust.door_ridges(edge_clust_all{i})));
                subs_clust_e.edge_ridges{i}  = unique(cell2mat(subs_clust.edge_ridges(edge_clust_all{i})));
                subs_clust_e.hall_ridges{i}  = unique(cell2mat(subs_clust.hall_ridges(edge_clust_all{i})));
                subs_clust_e.ridges2D{i}     = unique(cell2mat(subs_clust.ridges2D(edge_clust_all{i})));
                subs_clust_e.indices{i}      = edge_clust_inds{i};
                % Recompute edge indices
                [subs_clust_e.edgeind{i}, subs_clust_e.edgeEDT{i}] = ...
                               getEdgeData(subs_clust_e.indices{i}, voxels, ...
                                           shape, dx, subs_clust_e.center{i});
            else
                subs_clust_e.peaks.row{i}        = subs_clust.peaks.row{edge_clust_all{i}};
                subs_clust_e.peaks.ind{i}        = subs_clust.peaks.ind{edge_clust_all{i}};
                subs_clust_e.ridges1D{i}         = subs_clust.ridges1D{edge_clust_all{i}};
                subs_clust_e.ridges1D_lengths{i} = subs_clust.ridges1D_lengths{edge_clust_all{i}};
                subs_clust_e.skeleton{i}         = subs_clust.skeleton{edge_clust_all{i}};
                subs_clust_e.skeleton2D{i}       = subs_clust.skeleton2D{edge_clust_all{i}};
                subs_clust_e.center{i}           = subs_clust.center{edge_clust_all{i}};
                subs_clust_e.center_r1Dmin_diams{i} = subs_clust.center_r1Dmin_diams{edge_clust_all{i}};
                subs_clust_e.door_ridges{i}      = subs_clust.door_ridges{edge_clust_all{i}};
                subs_clust_e.edge_ridges{i}      = subs_clust.edge_ridges{edge_clust_all{i}};
                subs_clust_e.hall_ridges{i}      = subs_clust.hall_ridges{edge_clust_all{i}};
                subs_clust_e.ridges2D{i}         = subs_clust.ridges2D{edge_clust_all{i}};
                subs_clust_e.indices{i}          = edge_clust_inds{i};
                subs_clust_e.edgeind{i}          = subs_clust.edgeind{edge_clust_all{i}};
                subs_clust_e.edgeEDT{i}          = subs_clust.edgeEDT{edge_clust_all{i}};
            end
        end
        % Remove edge subs from subs_clust list
        rmv_edge = cell2mat(edge_clust_all);
        subs_clust.peaks.row(rmv_edge)        = [];
        subs_clust.peaks.ind(rmv_edge)        = [];
        subs_clust.ridges1D(rmv_edge)         = [];
        subs_clust.ridges1D_lengths(rmv_edge) = [];
        subs_clust.skeleton(rmv_edge)         = [];
        subs_clust.skeleton2D(rmv_edge)       = [];
        subs_clust.center(rmv_edge)           = [];
        subs_clust.center_r1Dmin_diams(rmv_edge) = [];
        subs_clust.door_ridges(rmv_edge)      = [];
        subs_clust.edge_ridges(rmv_edge)      = [];
        subs_clust.hall_ridges(rmv_edge)      = [];
        subs_clust.ridges2D(rmv_edge)         = [];
        subs_clust.indices(rmv_edge)          = [];
        subs_clust.edgeind(rmv_edge)          = [];
        subs_clust.edgeEDT(rmv_edge)          = [];

        % Update sub count
        numIntSubs = length(subs_clust.indices);
        num_subs   = numEdgeSubs + numIntSubs;
        edge_subs_log  = false(num_subs, 1);

        tElapsed = toc(tStart);
        % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Find edge of subunits:');

        time_log(timeLogIdx).Name = 'Find edge of subunits';
        time_log(timeLogIdx).Time = tElapsed;
        time_log(timeLogIdx).Units = 'sec';
        printTimeInfo(time_log(timeLogIdx));
        timeLogIdx = timeLogIdx + 1;

        %%%%%%%*****************************************************%%%%%%%
        %%%%%***************** DETERMINING NEIGHBORS *****************%%%%%
        %%%%%%%*****************************************************%%%%%%%
        tStart = tic;

        % Determine neighboring beads and touching bead pairs (crawl
        % spaces)
        neighbors_beads = cell(num_subs, 1);
        bead_pairs      = cell(num_subs, 1);
        % re-compute the subs_adjacencybeads matrix
        subs_adjacencybeads = cell(num_subs, 1);
        % Determine neighboring (adjacent) subunits using ridges2D
        subs_adjacencyr2D = cell(num_subs, 1);

        for h = 1 : num_subs
            % edge subs first, followed by interior subs
            if h <= numEdgeSubs
                subs_current = subs_clust_e;
                i = h;
            else
                subs_current = subs_clust;
                i = h - numEdgeSubs;
            end

            neighbors_beads{h} = unique(cell2mat(arrayfun(@(x) peaks.L7.beads{x}, ...
                                                                    subs_current.peaks.row{i}, ...
                                                                    'UniformOutput', false)));
            % creating 2-column array of: subunit number (col 1) & bead number (col 2)
            % for later adjacency matrix
            subs_adjacencybeads{h} = [h * ones(length(neighbors_beads{h}), 1), ...
                                      neighbors_beads{h}];

            % creating 2-column array of: subunit number (col 1) & ridge2D number (col 2)
            subs_adjacencyr2D{h} = [h * ones(length(subs_current.ridges2D{i}), 1), ...
                                    subs_current.ridges2D{i}];
        end

        % Sparse subunit-beads adjacency matrix
        subs_adjacencybeads = cell2mat(subs_adjacencybeads);
        subs_adjacencybeads = sparse(subs_adjacencybeads(:, 1), subs_adjacencybeads(:, 2), ...
                                     ones(size(subs_adjacencybeads, 1), 1), ...
                                     num_subs, num_beads);
        % Sparse subunit-ridges2D adjacency matrix
        subs_adjacencyr2D = cell2mat(subs_adjacencyr2D);
        subs_adjacencyr2D = sparse(subs_adjacencyr2D(:, 1), subs_adjacencyr2D(:, 2), ...
                                   ones(size(subs_adjacencyr2D, 1), 1), ...
                                   num_subs, length(ridges2D.indices));
        % Determine which subs share ridges2D by using subunit-r2D adjacency matrix
        neighbors_subsr2D = cell(num_subs, 1);
        crawl_spaces_log = true(num_r2D, 1); % defining crawl spaces between subunits
        for h = 1 : num_subs
            these_subs = [];
            these_r2D = find(subs_adjacencyr2D(h, :));
            for j = these_r2D(:)'
                adj_subs = find(subs_adjacencyr2D(:, j));
                these_subs = [these_subs; adj_subs];
                if isempty(adj_subs) || (numel(adj_subs) == 1 && adj_subs == h)
                    crawl_spaces_log(j) = false;
                end
            end
            these_subs = unique(these_subs);
            these_subs(these_subs == h) = [];
            % store subs
            neighbors_subsr2D{h} = these_subs;
        end
        ridges2D.crawl_log = crawl_spaces_log;

        % Determine neighboring (hallway) subunits using ridges1D
        neighbors_subsr1D = cell(num_subs, 1);
        % Create 3-column array for sparse matrix of subs graph
        subs_graph = cell(num_subs - 1, 1); % [sub #1, sub #2, door radius]
        for h = 1 : num_subs - 1
            % edge subs first, followed by interior subs
            if h <= numEdgeSubs
                subs_current_i = subs_clust_e;
                i = h;
            else
                subs_current_i = subs_clust;
                i = h - numEdgeSubs;
            end
            tracking = [];
            for hh = h + 1 : num_subs
                % edge subs first, followed by interior subs
                if hh <= numEdgeSubs
                    subs_current_j = subs_clust_e;
                    j = hh;
                else
                    subs_current_j = subs_clust;
                    j = hh - numEdgeSubs;
                end

                which_ridges = fastIntersect(subs_current_i.hall_ridges{i}, subs_current_j.hall_ridges{j}, 'elements');
                if ~isempty(which_ridges)
                    neighbors_subsr1D{h} = [neighbors_subsr1D{h}; hh];
                    neighbors_subsr1D{hh} = [neighbors_subsr1D{hh}; h];
                    for k = which_ridges(:)'
                        tracking = [tracking; [hh, ridges1D.doors{k}.radius]];
                    end
                end
            end
            neighbors_subsr1D{h} = unique(neighbors_subsr1D{h});
            subs_graph{h} = [h * ones(size(tracking, 1), 1), tracking];
        end
        subs_graph = subs_graph(~cellfun(@isempty, subs_graph));
        subs_graph = cell2mat(subs_graph);
        % Create sparse matrix representing graph of touching subunits, where
        % value in matrix is radius of connecting door
        if ~isempty(subs_graph)
            subs_adjacencysubs = sparse(subs_graph(:, 1), subs_graph(:, 2), ...
                                       subs_graph(:, 3), num_subs, num_subs);
        else
            subs_adjacencysubs = 0;
        end
        % Now create sparse adjacency matrix of subunits and doors (doors are
        % numbered as ridges1D are numbered)
        subsr1D_graph = [];
        for h = 1 : num_subs
            % edge subs first, followed by interior subs
            if h <= numEdgeSubs
                subs_current = subs_clust_e;
                i = h;
            else
                subs_current = subs_clust;
                i = h - numEdgeSubs;
            end
            % pair subunit with ridge number
            subsr1D_graph = [subsr1D_graph; ...
                             [i * ones(length(subs_current.door_ridges{i}), 1), ...
                              subs_current.door_ridges{i}]; ...
                             [i * ones(length(subs_current.edge_ridges{i}), 1), ...
                              subs_current.edge_ridges{i}]];
        end
        % add 3rd column of door radius
        subsr1D_graph = [sortrows(subsr1D_graph), zeros(size(subsr1D_graph, 1), 1)];
        for i = 1 : size(subsr1D_graph, 1)
            subsr1D_graph(i, 3) = ridges1D.doors{subsr1D_graph(i, 2)}.radius;
        end
        subs_adjacencydoors = sparse(subsr1D_graph(:, 1), subsr1D_graph(:, 2), ...
                                     subsr1D_graph(:, 3), num_subs, ridges1D.num);

        tElapsed = toc(tStart);
        % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Determine neighbors:');

        time_log(timeLogIdx).Name = 'Determine neighbors';
        time_log(timeLogIdx).Time = tElapsed;
        time_log(timeLogIdx).Units = 'sec';
        printTimeInfo(time_log(timeLogIdx));
        timeLogIdx = timeLogIdx + 1;

        %%%%%%%***************************************************%%%%%%%
        %%%%%**** LONGEST LENGTH OF SUBUNIT USING CONVEX HULL ***%%%%%
        %%%%%%%***************************************************%%%%%%%
        tStart = tic;
        convhull_longest = zeros(num_subs, 1);
        for h = 1 : num_subs
            % start with edge subs, followed by interior subs
            if h <= numEdgeSubs
                subs_current = subs_clust_e;
                i = h;
            else
                subs_current = subs_clust;
                i = h - numEdgeSubs;
            end
            try
                KK = convhull(voxels(subs_current.edgeind{i}, :));
                KK_inds = unique(KK);
                inds = subs_current.edgeind{i}(KK_inds);
                for ii = 1 : numel(inds) - 1
                    for jj = ii + 1 : numel(inds)
                        this_dist = norm(voxels(inds(ii),:) - voxels(inds(jj),:));
                        if this_dist > convhull_longest(h)
                            convhull_longest(h) = this_dist;
                        end
                    end
                end
            catch
                fprintf('%s\n', '??? Warning using ==> convhull');
                if numel(subs_current.edgeind{i}) <= 1
                    convhull_longest(h) = dx;
                else
                    inds = subs_current.edgeind{i};
                    for ii = 1 : numel(inds) - 1
                        for jj = ii + 1 : numel(inds)
                            this_dist = norm(voxels(inds(ii),:) - voxels(inds(jj),:));
                            if this_dist > convhull_longest(h)
                                convhull_longest(h) = this_dist;
                            end
                        end
                    end
                end
            end
        end
        tElapsed = toc(tStart);

        time_log(timeLogIdx).Name = 'Subunit length by convex hull';
        time_log(timeLogIdx).Time = tElapsed;
        time_log(timeLogIdx).Units = 'sec';
        printTimeInfo(time_log(timeLogIdx));
        timeLogIdx = timeLogIdx + 1;

        %%%%%%%****************************************************%%%%%%%
        %%%%%************* INTEGRIN-BINDING-PROTEIN MAP *************%%%%%
        %%%%%%%****************************************************%%%%%%%
        tStart = tic;
        % Use ligand map to later compute integrin-binding protein (RGD)
        % concentration within each subunit
        % Conversions:
        % 1 um^3 = 1e-15 L
        % 1 micromole = 6.02e17 molecules
        cell_diameter = 10;
        moving_window_diam = domain(2) / 60; % arbitrarily scaled to what we developed LOVAMAP on
        RGD_conc      = 500; % in micromoles/liter (per particle)
        % the RGD_cutoff value was chosen as 90% of the # of molecules per
        % voxel when dx = 2
        RGD_cutoff    = 2167200; % in # of molecules (old value:  1721136)
        % Generate ligand map, which gives avg # of RGD molecules at each voxel
        ligand_map  = ligandMap('shell', 'average', 'sphere', 'all', bead_struct, ...
                                moving_window_diam, RGD_conc, vsVoxel_log, shape, dx);

        % TODO Janky crop version to temporarily produce better hotspot data
        ligand_map_crop = ligand_map;
        ligand_map_crop(:, :, (size(ligand_map, 3) - (shell_thickness + 1)*dx) : end) = 0;

        % Used for figuring out what voxels remain in the convex hull
        cropHeight = domain(6) - dx * shell_thickness;

        % Figure out what voxels still remain in the convex hull after crop
        cHullXYZ = voxels(vsVoxel_log, :);
        cHullXYZ_crop = cHullXYZ(cHullXYZ(:, 3) < cropHeight, :);
        nVoxelsInCroppedHull = size(cHullXYZ_crop, 1);

        % Figure out what bead voxels still remain after the crop
        allBeadsXYZ = voxels(data.allBeads, :);
        allBeadsXYZ_crop = allBeadsXYZ(allBeadsXYZ(:, 3) < cropHeight, :);
        nVoxelsInCroppedBeads = size(allBeadsXYZ_crop, 1);

        hotspot_inds_crop = ligand_map_crop >= RGD_cutoff;
        hotspot_total_crop = sum(ligand_map_crop(hotspot_inds_crop), 'all') / 6.02e17;
        hotspot_ratio_crop = sum(hotspot_inds_crop, 'all') / (nVoxelsInCroppedHull + nVoxelsInCroppedBeads);

        % Original way of computing hotspot stuff
        % hotspot_inds  = ligand_map >= RGD_cutoff;
        % hotspot_total = sum(ligand_map(hotspot_inds) 'all') / 6.02e17; % in micromoles
        % % for ratio, divide by all scaffold voxels (i.e., void + beads)
        % hotspot_ratio = sum(hotspot_inds, 'all') / (num_VsVoxels + numel(bead_struct.AllBeads));

        hotspot_ratio = hotspot_ratio_crop;
        data.ligandMap = ligand_map;
        data.LigandHotspots = hotspot_inds_crop;

        % Isolate ligand map of surface subunits
        for i = 1 : num_2Dedgesubs
            tot_edgeLigand         = sum(ligand_map(edge_2Dsubs{i}.indices)); % # molecules
            edgesub_SA             = round2sigdig(length(edge_2Dsubs{i}.indices) * dx^2, dx); % µm^2
            edge_2Dsubs{i}.RGDconc = tot_edgeLigand / (edgesub_SA * 6.02e17);
        end

        tElapsed = toc(tStart);
        % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Ligand concentration:');

        time_log(timeLogIdx).Name = 'Ligand concentration';
        time_log(timeLogIdx).Time = tElapsed;
        time_log(timeLogIdx).Units = 'sec';
        printTimeInfo(time_log(timeLogIdx));
        timeLogIdx = timeLogIdx + 1;

        %%%%%%%***************************************************%%%%%%%
        %%%%%***************** SUBUNIT INFORMATION *****************%%%%%
        %%%%%%%***************************************************%%%%%%%
        % Create a cell array of subunit structs
        tStart = tic;
        subunits = cell(num_subs, 1);
        max_numpksE = 0;
        % NOTE!!: Remember i tracks the subunit number, but j tracks the
        % index into the data used to fill the subunits!!!

        for i = 1 : num_subs
            % edge subs first, followed by interior subs
            if i <= numEdgeSubs
                subs_final = subs_clust_e;
                j = i;
            else
                subs_final = subs_clust;
                j = i - numEdgeSubs;
            end
            subunits{i} = struct;
            % subunit #
            subunits{i}.num = i;
            % unique identifier
            subunits{i}.uniqueID = [num2str(subunits{i}.num) '-' generateRandomID()];
            % all voxel indices that make up subunit
            subunits{i}.indices = subs_final.indices{j};
            % edge voxel indices
            subunits{i}.edgeindices = subs_final.edgeind{j};
            % for plotting
            subunits{i}.gridpts = voxels(subunits{i}.indices, :);
            % center coordinate
            subunits{i}.centerCoord = mean(subunits{i}.gridpts, 1);
            % voxel indices of peaks
            subunits{i}.peaks = subs_final.peaks.ind{j};
            % EDT at peaks
            subunits{i}.EDTofPeaks = peaks.L7.EDT(subs_final.peaks.row{j});
            % largest enclosed sphere (diameter)
            subunits{i}.largestSphereDiam = max(subunits{i}.EDTofPeaks) * 2;
            % largest enclosed sphere (volume in pL)
            subunits{i}.largestSphereVol = (4/3) * pi * (subunits{i}.largestSphereDiam / 2).^3;
            % number of peaks
            subunits{i}.numPeaks = length(subunits{i}.peaks);
            % ridge1D #s
            subunits{i}.ridges1D = subs_final.ridges1D{j};
            % ridge1D lengths per ridge
            subunits{i}.ridge1DLengths = subs_final.ridges1D_lengths{j};
            % ridge2D #s
            subunits{i}.ridges2D = subs_final.ridges2D{j};
            % skeleton voxels
            subunits{i}.skeletonVoxs = subs_final.skeleton{j};
%             % create alphaShape from points that make up subunit
%             subunits{i}.alphShape = alphaShape(subunits{i}.gridpts);
            subunits{i}.skeleton2DVoxs = subs_final.skeleton2D{j};
            subunits{i}.centerSkeleton = subs_final.center{j};
            % subunit door data
            subunits{i}.doorRidges = subs_final.door_ridges{j};
            % subunit edge peaks
%             subunits{i}.edgeRidges = subs_final.edge_ridges{j};
            if i <= numEdgeSubs
                subunits{i}.edgePeaksInd = edge_clust_pks{j};
                % surface area
                outer_edge = fastIntersect(subunits{i}.edgeindices, edge_ind, 'elements');
                subunits{i}.outerSurfInds = outer_edge;
                subunits{i}.outerSurfArea = round2sigdig(length(outer_edge) * dx^2, dx); % in um^2
            else
                subunits{i}.edgePeaksInd  = [];
                subunits{i}.outerSurfInds = [];
                subunits{i}.outerSurfArea = [];
            end
            % subunit volume
            subunits{i}.volume = round2sigdig(length(subunits{i}.indices) * dx^3, dx) / 1000; % in pL
            % surface area
            subunits{i}.surfArea = round2sigdig(length(subunits{i}.edgeindices) * dx^2, dx) / 1000; % in um^2/1000
            % characteristic length
            if subunits{i}.surfArea < dx * 5e-05 % avoid dividing by 0
                subunits{i}.charLength = 0;
            else
                subunits{i}.charLength = subunits{i}.volume / subunits{i}.surfArea;
            end
            % end-to-end length using convex hull
            subunits{i}.convLength = convhull_longest(i);
            % average 'interior' door diameter of subunit (um)
            if numel(subunits{i}.EDTofPeaks) == 1
                subunits{i}.avgInternalDiam = subunits{i}.EDTofPeaks * 2;
            else
                % subunits{i}.narrowestDiam = mean(subs_final.center_r1Dmin_diams{j});
                subunits{i}.avgInternalDiam = mean(data.EDT(subunits{i}.centerSkeleton) * 2);
            end
            % bead neighbors
            subunits{i}.beadNeighbors = neighbors_beads{i};
            % subunit neighbors that are connected by hallways (ridges1D)
            subunits{i}.subNeighborsr1D = neighbors_subsr1D{i};
            % subunit neighbors that are adjacent (ridges2D)
            subunits{i}.subNeighborsr2D = neighbors_subsr2D{i};
            % # of hallways leaving subunit
            subunits{i}.numHalls = numel(subs_final.hall_ridges{j}) + numel(subs_final.edge_ridges{j});
            % # of crawl spaces
            subunits{i}.numCrawlSpaces = numel(subunits{i}.ridges2D(ridges2D.crawl_log(subunits{i}.ridges2D)));
            % # of connected pores
            subunits{i}.numConnectPores = numel(subunits{i}.subNeighborsr1D);
            % # of surrounding pores
            subunits{i}.numSurroundPores = numel(subunits{i}.subNeighborsr2D);
            % ratio of halls to crawl spaces (normalized connectivity)
            subunits{i}.normNeigh = subunits{i}.numHalls / subunits{i}.numCrawlSpaces;
            % average of door radii
            subunits{i}.avgDoorDiam = mean(arrayfun(@(x) ridges1D.doors{x}.radius * 2, subunits{i}.doorRidges));
            subunits{i}.avgDoorDiam(isnan(subunits{i}.avgDoorDiam)) = 0;
            % largest door diameter
            subunits{i}.largestDoorDiam = max(arrayfun(@(x) ridges1D.doors{x}.radius * 2, subunits{i}.doorRidges));
            subunits{i}.largestDoorDiam(isempty(subunits{i}.largestDoorDiam)) = 0;
            % largest door sphere volume
            subunits{i}.largestDoorVol = (4/3) * pi * (subunits{i}.largestDoorDiam / 2).^3;
            subunits{i}.largestDoorVol(isempty(subunits{i}.largestDoorVol)) = 0;
            % smallest door diameter
            subunits{i}.smallestDoorDiam = min(arrayfun(@(x) ridges1D.doors{x}.radius * 2, subunits{i}.doorRidges));
            subunits{i}.smallestDoorDiam(isempty(subunits{i}.smallestDoorDiam)) = 0;
            % smallest door sphere volume
            subunits{i}.smallestDoorVol = (4/3) * pi * (subunits{i}.smallestDoorDiam / 2).^3;
            subunits{i}.smallestDoorVol(isempty(subunits{i}.smallestDoorVol)) = 0;
            % edge subunit or not
%             if fastIntersect(subunits{i}.edgeindices, edge_ind, 'bool')
%                 subunits{i}.edge = true;
%                 edge_subs_log(i) = true;
%             else
%                 subunits{i}.edge = false;
%             end
            if i <= numEdgeSubs
                subunits{i}.edge = true;
                edge_subs_log(i) = true;
                % max number of peaks of edge subunits
                max_numpksE = max(max_numpksE, numel(subunits{i}.peaks));
            else
                subunits{i}.edge = false;
                edge_subs_log(i) = false;
            end

            %%%%% MEAN LOCAL THICKNESS %%%%%
            % Compute mean local thickness
            if combine_edge_subs && subunits{i}.edge == true
                subunits{i}.meanLocalThickness = NaN;
            else
                if subunits{i}.volume * dx > 400
                    subunits{i}.meanLocalThickness = NaN;
                else
                    local_thick_mat = localThickness(subunits{i}.indices, ...
                                      subs_final.skeleton2D{j}, data.EDT(subs_final.skeleton2D{j}), ...
                                      voxels, dx, shape);
                    subunits{i}.meanLocalThickness = (sum(local_thick_mat) * dx^3) / (subunits{i}.volume * 1000);
                end
            end

            %%%%%%% EDGE EDT AND INTEGRIN STUFF %%%%%%%
            % EDT (from center) along edge of subunit
            subunits{i}.edgeEDT = subs_final.edgeEDT{j};
            % total integrin-binding protein concentration within subunit
            % using conversions:
            % 1 um^3 = 1e-15 L
            % 1 micromole = 6.02e17 molecules
            tot_ligand               = sum(ligand_map(subunits{i}.indices)); % # molecules
            sub_volume               = subunits{i}.volume * 1000; % liters
            subunits{i}.RGDconc      = (tot_ligand * 1e15) / (sub_volume * 6.02e17);

            %%%%%% 'PARTICLE EDGES' OF A PORE %%%%%%
            % Descriptor for pore 'edges' - i.e., number of surrounding
            % particle + 1D-ridge pairs per pore
            total_edges = 0;
            for k = subunits{i}.doorRidges(:)'
                %if data.EDT(ridges1D.min(k)) > (dx*2) % check if the particle is substantially far from the 1D-ridge (otherwise there's really no 'edge' to consider..)
                these_pairs = fastIntersect(ridges1D.beads{k}, subunits{i}.beadNeighbors, 'size');
                % For each 1D-ridges exiting the pore, we're counting
                % up the number of particles surrounding the 1D-ridges
                % that are also surrounding the pore
                total_edges = total_edges + these_pairs;
                %end
            end
            subunits{i}.numFauxEdges = total_edges;
        end

        tElapsed = toc(tStart);
        % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Accumulate subunit data:');

        time_log(timeLogIdx).Name = 'Second set of pore data, inc. mean local thickness';
        time_log(timeLogIdx).Time = tElapsed;
        time_log(timeLogIdx).Units = 'sec';
        printTimeInfo(time_log(timeLogIdx));
        timeLogIdx = timeLogIdx + 1;

        tStart = tic;
        %%%%% VVF OF INTERIOR PORES ONLY %%%%%
        % Compute void volume fraction of interior pores only
        int_inds = 0;
        for i = 1 : num_subs
            if ~subunits{i}.edge
                int_inds = int_inds + numel(subunits{i}.indices);
            end
        end
        porosity_int = int_inds / sum(in_log);

        %%%%%% PATH STUFF %%%%%%%
        tic;
        % Compute pores per path by locating pores that share peaks along path(s)
        data.paths.subunits = cell(data.paths.num_paths, 1);
        data.paths.subunits_key = cell(data.paths.num_paths, 1); % 'interior' pores is 'true', 'exterior' pores is 'false'
        data.paths.surrounding_particles.int = cell(data.paths.num_paths, 1);
        data.paths.surrounding_particles.ext = cell(data.paths.num_paths, 1);
        % Generate matrix of path pore/door data:
        % Rows refer to paths
        % Columns list alternating: pore data, door diameter, pore data, door diameter,
        % etc. going from interior of scaffold outward along path
        data.paths.path_matrix_pores = cell(data.paths.num_paths, 1);
        data.paths.path_matrix_nodes = cell(data.paths.num_paths, 1);
        % Find first pore from center node
        for i = 1 : numel(subunits)
            if fastIntersect(peaks.L7.indices(data.paths.center_pk), subunits{i}.peaks, 'bool')
                center_pore = i;
                break;
            end
        end
        pore_list = 1 : numel(subunits);

        % Find remaining pores from 1D-ridges
        for i = 1 : data.paths.num_paths
            remaining_pores = true(1, numel(subunits));
            row_cntr = 1; % for path_matrix
            % Grab r1D data for the path
            these_r1Ds = data.paths.path_r1Ds{i};

            % First pore will always be the center pore for each path
            data.paths.subunits{i} = center_pore;
            data.paths.subunits_key{i} = ~subunits{center_pore}.edge;
            if subunits{center_pore}.edge
                    % exterior pores
                    data.paths.surrounding_particles.ext{i} = subunits{center_pore}.beadNeighbors;
                else
                    % interior pores
                    data.paths.surrounding_particles.int{i} = subunits{center_pore}.beadNeighbors;
            end
            % fill matrix cell with pore data
            data.paths.path_matrix_pores{i}{row_cntr} = [subunits{center_pore}.convLength, subunits{center_pore}.avgInternalDiam, subunits{center_pore}.largestSphereDiam];
            row_cntr = row_cntr + 1;

            % Find remaining pores and doors in order
            last_pore = center_pore;
            remaining_pores(last_pore) = false;
            for j = these_r1Ds(:)'
                this_pore = [];
                for k = pore_list(remaining_pores)
                    if fastIntersect(j, subunits{k}.ridges1D, 'bool')
                        % since we removed the other pore (connected to the ridge) from the list of pores, the ridge
                        % should only be associated with one other pore max
                        this_pore = k;
                        break;
                    end
                end
                if ~isempty(this_pore)
                    % Fill path matrix with door data
                    data.paths.path_matrix_pores{i}{row_cntr} = ridges1D.doors{j}.radius * 2;
                    row_cntr = row_cntr + 1;

                    % Fill in pore and door info
                    data.paths.subunits{i} = [data.paths.subunits{i}; this_pore];
                    data.paths.subunits_key{i} = [data.paths.subunits_key{i}; ~subunits{this_pore}.edge];
                    if subunits{this_pore}.edge
                            % exterior pores
                            data.paths.surrounding_particles.ext{i} = [data.paths.surrounding_particles.ext{i}; subunits{this_pore}.beadNeighbors];
                        else
                            % interior pores
                            data.paths.surrounding_particles.int{i} = [data.paths.surrounding_particles.int{i}; subunits{this_pore}.beadNeighbors];
                    end
                    % fill matrix cell with pore data
                    data.paths.path_matrix_pores{i}{row_cntr} = [subunits{this_pore}.convLength, subunits{this_pore}.avgInternalDiam, subunits{this_pore}.largestSphereDiam];
                    row_cntr = row_cntr + 1;

                    % Update last pore to remove it from potential-pores list. Add the previous pore back to the list to
                    % accommodate large pores that weave in and out of the path
                    remaining_pores(last_pore) = true;
                    remaining_pores(this_pore) = false;
                    last_pore = this_pore;
                end
            end
            data.paths.surrounding_particles.ext{i} = unique(data.paths.surrounding_particles.ext{i});
            data.paths.surrounding_particles.int{i} = unique(data.paths.surrounding_particles.int{i});

            % Alternative matrix output of nodes and doors
            % % Check for length discrepancy
            % if numel(data.paths.path_r1Ds{i}) > numel(data.paths.node_pks{i}) || ...
            %   numel(data.paths.node_pks{i}) > numel(data.paths.path_r1Ds{i}) + 1
            %     error('The length of path_r1Ds is greater than the length of node_pks. There is an error with the data.');
            % end

            % Determine length of shorter vector
            minLength = numel(data.paths.path_r1Ds{i});
            path_matrix_row = zeros(1, 2 * minLength);

            % Interleave elements from node_pks and path_r1Ds
            path_matrix_row(1:2:end) = data.EDT(peaks.L7.indices(data.paths.node_pks{i}(1 : minLength))) * 2;
            path_matrix_row(2:2:end) = arrayfun(@(x) ridges1D.doors{x}.radius * 2, data.paths.path_r1Ds{i}(1 : minLength));

            % Append any remaining elements from node_pks
            if length(data.paths.node_pks{i}) > minLength
                path_matrix_row = [path_matrix_row, data.paths.node_pks{i}(minLength + 1 : end)]; % shouldn't be more than a single value..
            end

            data.paths.path_matrix_nodes{i} = path_matrix_row;
        end

        tElapsed = toc(tStart);

        time_log(timeLogIdx).Name = 'Third set of pore data, inc. path data';
        time_log(timeLogIdx).Time = tElapsed;
        time_log(timeLogIdx).Units = 'sec';
        printTimeInfo(time_log(timeLogIdx));
        timeLogIdx = timeLogIdx + 1;

        %
        %
        %
        %
        %
        %     num_nodes = numel(data.paths.node_pks{i});
        %     these_path_subs = zeros(num_nodes, 1);
        %     boolean_key = zeros(num_nodes, 1);
        %     % we don't know ahead of time how many surrounding particles
        %     these_surr_particles_int = [];
        %     these_surr_particles_ext = [];
        %     % Order of how we search unfortunately matters here because we
        %     % need to keep track of the order of pores along path
        %     for j = 1 : num_nodes
        %         for k = 1 : numel(subunits)
        %             if fastIntersect(peaks.L7.indices(data.paths.node_pks{i}(j)), ...
        %                                     subunits{k}.peaks, 'bool')
        %                 these_path_subs(j) = k;
        %                 if subunits{k}.edge
        %                     % exterior pores
        %                     boolean_key(j) = false;
        %                     these_surr_particles_ext = [these_surr_particles_ext; subunits{j}.beadNeighbors];
        %                 else
        %                     % interior pores
        %                     boolean_key(j) = true;
        %                     these_surr_particles_int = [these_surr_particles_int; subunits{j}.beadNeighbors];
        %                 end
        %                 break;
        %             end
        %         end
        %     end
        %     data.paths.subunits{i} = unique(these_path_subs, 'stable');
        %     data.paths.subunits_key{i} = boolean_key;
        %     data.paths.surrounding_particles.int{i} = unique(these_surr_particles_int);
        %     data.paths.surrounding_particles.ext{i} = unique(these_surr_particles_ext);
        % end
        %
        % % Generate matrix of path pore/door data:
        % % Rows refer to paths
        % % Columns list alternating: pore data, door diameter, pore data, door diameter,
        % % etc. going from interior of scaffold outward along path
        % path_matrix = cell(data.paths.num_paths, 1);
        % for i = 1 : data.paths.num_paths
        %     % remove extra cells after
        %     num_pores = numel(data.paths.subunits{i});
        %     path_r1Ds = sort(data.paths.path_r1Ds{i});
        %     path_matrix_row = cell(num_pores * 3, 1);
        %     end_now = 0;
        %     cntr = 1;
        %     large_pore_r1D_trkr = data.paths.path_r1Ds{i};
        %     for j = 1 : num_pores
        %         if end_now == 1 % theoretically this check isn't needed
        %             break;
        %         end
        %         current_pore = data.paths.subunits{i}(j);
        %         % fill matrix cell with pore data
        %         path_matrix_row{cntr} = [subunits{current_pore}.convLength, subunits{i}.avgInternalDiam, subunits{i}.largestSphereDiam];
        %         % fill next matrix cell with door diameter
        %         cntr = cntr + 1;
        %         if j < num_pores
        %             current_pore_doors = sort(subunits{current_pore}.doorRidges);
        %             next_pore_doors = sort(subunits{data.paths.subunits{i}(j+1)}.doorRidges);
        %             these_doors = fastIntersect(current_pore_doors, next_pore_doors, 'elements');
        %             this_door = fastIntersect(sort(these_doors), path_r1Ds, 'elements');
        %             if numel(this_door) > 1
        %                 % Examples where path goes in and out of one larger pore
        %                 [drs, inds, ~] = intersect(large_pore_r1D_trkr, this_door, 'stable');
        %                 this_door = drs(1);
        %                 large_pore_r1D_trkr(inds(1)) = [];
        %             end
        %             path_matrix_row{cntr} = ridges1D.doors{this_door}.radius * 2;
        %         else
        %             this_door = path_r1Ds(end);
        %             path_matrix_row{cntr} = ridges1D.doors{this_door}.radius * 2;
        %             end_now = 1;
        %         end
        %         cntr = cntr + 1;
        %     end
        %     path_matrix_row = path_matrix_row(~cellfun(@isempty, path_matrix_row));
        %     path_matrix{i} = path_matrix_row;
        % end
        % data.paths.path_matrix = path_matrix;

        tStart = tic;
        %%%%% ELLIPSOIDS %%%%%
        % Get ellipsoid/orientation/isotropy/directionality data
        % Column 1 : 3 of ellipse_data gives the eigenvalues of the PCA for
        % the subunit point cloud. The radii of the corresponding
        % ellipsoid are chosen as the squareroot of the eigenvalues (= the SV), meaning
        % the lengths of the ellipsoid axes are sqrt(eg's) * 2
        ellipse_data = ellipseData(subunits, 0, 0, 0, 0);
        ellipse_class = zeros(num_subs, 1);
        for i = 1 : num_subs
            % convert eigenvalues to ellipsoid lengths
            subunits{i}.ellipse_lengths = sqrt(ellipse_data(i, 1 : 3)) .* 2;
            % By definition, A >= B >= C
            A = subunits{i}.ellipse_lengths(1);
            B = subunits{i}.ellipse_lengths(2);
            C = subunits{i}.ellipse_lengths(3);
            if sum(subunits{i}.ellipse_lengths < dx * 5e-05) > 0 % problems with nearly-0 values
                subunits{i}.ellipse_lengths = zeros(1, 3);
                A = 0;
                B = 0;
                C = 0;
            end
            % compute ellipsoid-shape metric

            % KEY: [0, 1, 2, 3] => [undefined, sphere, pancake, tube]

            if sum(subunits{i}.ellipse_lengths == 0)
                % CASE 0: Empty (set to 0)
                ellipse_class(i) = 0;
            elseif C / A >= 0.8
                % Case 1: SPHERE = 1
                % axis 3 is at least 80% the length of axis 1 (i.e., all axes are similar in length)
                ellipse_class(i) = 1;
            elseif (B - C) / (A - C) > 0.5 && C / A < (1/3)
                % Case 2: PANCAKE = 2
                % axis 1 and 2 are similar relative to axis 3, and
                % axis 3 is under 1/3 the length of axis 1
                ellipse_class(i) = 2;
            elseif (B - C) / (A - C) < 0.5 && B / A < (1/3)
                % Case 3: TUBE = 3
                % axis 2 and 3 are similar relative to axis 1, and
                % axis 2 is under 1/3 the length of axis 1
                ellipse_class(i) = 3;
            else
                % Case 4: Undefined (set to 0)
                ellipse_class(i) = 0;
            end
            subunits{i}.ellipse_class = ellipse_class(i);
            % isotropy
            subunits{i}.isotropy = abs(ellipse_data(i, 4));
        end

        %%%%%% MORE LIGAND STUFF %%%%%%%
        % Compute accessible ligand (in micromoles) by identifying edge
        % voxels of neighbor beads and multiplying by thickness of shell
        % Create 3D matrix where beads are identified by number
        beads_by_number = zeros(shape);
        for i = 1 : num_beads
            beads_by_number(data.beads{i}) = i;
        end
        beads_by_number_padded = padarray(beads_by_number, [1 1 1], 0, 'both');
        for i = 1 : num_subs
            num_edgevoxs = numel(subunits{i}.edgeindices);
            edge_shellbool = false(num_edgevoxs, 1);
%             edge_beadnum   = zeros(num_edgevoxs, 1);
            edge_shellinds = zeros(num_edgevoxs, 1);
            for j = 1 : num_edgevoxs
                bounce = false;
                % check if edge voxel neighbors a bead
                [x,y,z] = ind2sub(shape, subunits{i}.edgeindices(j));
                % shift indices to adjust for pad
                x = x + 1;
                y = y + 1;
                z = z + 1;
                for k = x - 1 : x + 1
                    for l = y - 1 : y + 1
                        for m = z - 1 : z + 1
                            if beads_by_number_padded(k,l,m) ~= 0
                                edge_shellbool(j) = true;
%                                 edge_beadnum(j)   = beads_by_number(k,l,m);
                                edge_shellinds(j) = sub2ind(shape, k - 1, l - 1, m - 1);
                                bounce = true;
                            end
                            if bounce
                                break;
                            end
                        end
                        if bounce
                            break;
                        end
                    end
                    if bounce
                        break;
                    end
                end
            end
            edge_shellinds = edge_shellinds(edge_shellinds ~= 0);
            subunits{i}.touchingBeadindices = edge_shellinds;
            % Store edge ligand information
            num_edgeshellvoxs = sum(edge_shellbool);
            num_edgeopenvoxs = sum(~edge_shellbool);
            SA_edgeshell = round2sigdig(num_edgeshellvoxs * dx^2, dx);
            SA_edgeopen  = round2sigdig(num_edgeopenvoxs * dx^2, dx);
            tot_equiv_vol = SA_edgeshell * shell_thickness; % in um^3
            % using conversion:
            % 1 um^3 = 1e-15 L
            subunits{i}.RGDaccessible = RGD_conc * tot_equiv_vol * 1e-15;
            % Storing particle surface area of pore and open (non-particle) surface area of pore
            subunits{i}.SA_beads = SA_edgeshell;
            subunits{i}.SA_open  = SA_edgeopen;
        end


        %%%%%% EDGE SUBUNIT STUFF - TO-DO %%%%%%
        % Compute 'edge subunit' descriptors

        % [To-do]


        %%%%% EIGNEVALUES %%%%%
        % Get largest eigenvalue of global adjacency matrices
        % beads
        bead_adj_ones = beads_adjacency;
        bead_adj_ones(bead_adj_ones ~= 0) = 1;
        [beads_beads_maxe, ~] = power_method(bead_adj_ones + bead_adj_ones', 5000, 1e-12);
        % peaks
        peak_adj_ones = peaks.adjacency;
        peak_adj_ones(peak_adj_ones ~= 0) = 1;
        [peaks_peaks_maxe, ~] = power_method(peak_adj_ones + peak_adj_ones', 5000, 1e-12);
        % beads
        subs_adj_ones = subs_adjacencysubs;
        subs_adj_ones(subs_adj_ones ~= 0) = 1;
        [subs_subs_maxe, ~]   = power_method(subs_adj_ones + subs_adj_ones', 5000, 1e-12);
%         subs_beads  = subs_adjacencybeads;
%         subs_r2D    = subs_adjacencyr2D;
%         subs_doors  = subs_adjacencydoors;

        %%%%% INTERIOR PORES PRE SURROUNDING PARTICLES %%%%%%
        % Compute # interior subs / # particles surrounding interior subs
        interior_beads = unique(cell2mat(arrayfun(@(x) subunits{x}.beadNeighbors, ...
                                find(~edge_subs_log)', 'UniformOutput', false)'));
        subs2beads = (num_subs - sum(edge_subs_log)) / numel(interior_beads);

        tElapsed = toc(tStart);
        % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Compute eigenvalues of adjacency:');

        time_log(timeLogIdx).Name = 'Third set of pore data';
        time_log(timeLogIdx).Time = tElapsed;
        time_log(timeLogIdx).Units = 'sec';
        printTimeInfo(time_log(timeLogIdx));
        timeLogIdx = timeLogIdx + 1;

        % Gather descriptor data
        tStart = tic;

        % Gather doors data
        interior_door_ridges = unique(cell2mat(arrayfun(@(x) subunits{x}.doorRidges, ...
                                               1 : num_subs, 'UniformOutput', false)'));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%******************** DESCRIPTORS *********************%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store data by descriptors
        data.Descriptors           = struct;
        data.Descriptors.Global    = struct;
        data.Descriptors.NonSubs = struct;
        data.Descriptors.Subs      = struct;
        % Fill global descriptors
        data.Descriptors.Global.names =           {'dx';
                                                   'voxels';
                                                   '# Particles';
                                                   '# Particle Contacts';
                                                   'Void Vol Fraction';
                                                   'Particle Fraction';
                                                   'Void Area Fraction';
                                                   'Void Vol Fraction of Interior Pores'
                                                   '# 3D Pores';
                                                   '# Interior 3D Pores';
                                                   '# Int 3D-Pores / # Particles Surr Int Pores';
                                                   '# Exterior Doors';
                                                   '# Interior Doors';
                                                   '# Paths';
                                                   'Particle Adj Max Eig';
                                                   'Peak Adj Max Eig';
                                                   'Pore Adj Max Eig';
                                                   'Max # Pks of Edge Pore';
                                                   'Ligand Hotspots Volume Fraction';
                                                   'Max # Equidistant Particles'};

        data.Descriptors.NonSubs.names =         {'Particle Diameter (um)';
                                                   'Particle Coordination #';
                                                   'Touching Particle Coord #';
                                                   'Exit Door Diameter (um)';
                                                   'Internal Door Diameter (um)';
                                                   'Crawl Space Width (um)';
                                                   'Path Length (um)';
                                                   'Tortuosity by Length';
                                                   'Tortuosity by Volume (pL)';
                                                   '# Interior 3D Pores Traversed by Path';
                                                   '# Exterior 3D Pores Traversed by Path';
                                                   '# Particles Enclosing Interior Path';
                                                   '# Particles Enclosing Exiting Path';
                                                   '# Bottlenecks (total)';
                                                   'Avg Bottleneck Diameter (total) (um)';
                                                   '# Bottlenecks (doors)';
                                                   'Avg Bottleneck Diameter (doors) (um)';
                                                   'Surface Ligand Conc (µmoles / µm^2)';
                                                   'Surface Accessible Ligand (µmoles)';
                                                   'Region Vols, <1 um (pL)';
                                                   'Region Vols, 10 um (pL)';
                                                   'Region Vols, 30 um (pL)';
                                                   'Region Vols, 60 um (pL)'};

        data.Descriptors.Subs.names =             {'Volume (pL)';
                                                   'Surface Area (um2 / 1000)';
                                                   'Characteristic Length (um)';
                                                   'Longest Length (um)';
                                                   'Avg Internal Diam (um)';
                                                   'Aspect Ratio';
                                                   '# Hallways';
                                                   '# Crawl Spaces';
                                                   '# Connected Pores';
                                                   '# Surrounding Pores';
                                                   '# Particle Edges'
                                                   'Normalized Neighbors';
                                                   'Mean Local Thickness (um)';
                                                   '# Peaks';
                                                   '# Surrounding Particles';
                                                   'Largest Enclosed Sphere Diameter (um)';
                                                   'Largest Door Diameter (um)';
                                                   'Smallest Door Diameter (um)';
                                                   'Ellipsoid Axis 1 Length (um)';
                                                   'Ellipsoid Axis 2 Length (um)';
                                                   'Ellipsoid Axis 3 Length (um)';
                                                   'Ellipsoid Classification';
                                                   'Isotropy';
                                                   'Ligand Concentration (umoles / L)';
                                                   'Accessible Ligand (umoles)';
                                                   'Particle Surface Area (um^2)';
                                                   'Open (Non-Particle) Surface Area (um^2)';
                                                   'x Centroid';
                                                   'y Centroid';
                                                   'z Centroid';
                                                   'UniqueID'};

        % Fill global descriptors
        data.Descriptors.Global.dx                  = dx;
        data.Descriptors.Global.voxels              = nVoxels;
        data.Descriptors.Global.numBeads            = num_beads;
        data.Descriptors.Global.numContacts         = sum(touching_beads_log);
        data.Descriptors.Global.voidVolFract        = void_vol_fract;
        data.Descriptors.Global.particleVolFract    = particle_fract;
        data.Descriptors.Global.voidAreaFract       = void_area_fract;
        data.Descriptors.Global.voidVolFractInt     = porosity_int;
        data.Descriptors.Global.numSubs             = num_subs;
        data.Descriptors.Global.numIntSubs          = num_subs - sum(edge_subs_log);
        data.Descriptors.Global.subs2Beads          = subs2beads;
        data.Descriptors.Global.numDoors_exterior   = peaks_e.num;
        data.Descriptors.Global.numDoors_interior   = length(interior_door_ridges);
        data.Descriptors.Global.numPaths            = data.paths.num_paths;
        data.Descriptors.Global.beadEigenvalue      = beads_beads_maxe;
        data.Descriptors.Global.peakEigenvalue      = peaks_peaks_maxe;
        data.Descriptors.Global.subEigenvalue       = subs_subs_maxe;
        data.Descriptors.Global.maxNumPksEdge       = max_numpksE;
        data.Descriptors.Global.RGDhotspotsRatio    = hotspot_ratio;
        data.Descriptors.Global.maxEquidistBeads    = max(peaks.L7.numBeads);
        % Fill non-subunit descriptors
        % >> beads
        data.Descriptors.NonSubs.beadDiams          = bead_diams;
        data.Descriptors.NonSubs.beadCoord          = bead_struct.CoordCount;
        data.Descriptors.NonSubs.beadTouch          = bead_struct.TouchingCount;
        % >> doors & crawl spaces
        data.Descriptors.NonSubs.diamDoors_exterior = cellfun(@(x) x * 2, peaks_e.radius);
        data.Descriptors.NonSubs.diamDoors_interior = sort(arrayfun(@(x) ridges1D.doors{x}.radius * 2, ...
                                                        interior_door_ridges));
        data.Descriptors.NonSubs.crawlWidth         = ridges2D.crawl_widths(crawl_spaces_log);
        % >> paths
        data.Descriptors.NonSubs.pathLengths        = data.paths.path_lengths;
        data.Descriptors.NonSubs.pathTortuosity_lin = data.paths.tortuosity.linear;
        data.Descriptors.NonSubs.pathTortuosity_vol = data.paths.tortuosity.volume / 1000;
        data.Descriptors.NonSubs.pathNumPores_int = cellfun(@(x, y) numel(x(y)), data.paths.subunits, data.paths.subunits_key);
        data.Descriptors.NonSubs.pathNumPores_ext   = cellfun(@(x, y) numel(x(~y)), data.paths.subunits, data.paths.subunits_key);
        data.Descriptors.NonSubs.pathNumSurrBeads_int = cellfun(@(x) numel(x), data.paths.surrounding_particles.int);
        data.Descriptors.NonSubs.pathNumSurrBeads_ext = cellfun(@(x) numel(x), data.paths.surrounding_particles.ext);
        data.Descriptors.NonSubs.pathNumBotnecksTot   = cellfun(@(x) numel(x), data.paths.neck_data);
        data.Descriptors.NonSubs.pathAvgBothecksTot   = cellfun(@(x) mean(x), data.paths.neck_data);
        data.Descriptors.NonSubs.pathNumBotnecksDoor  = cellfun(@(x) numel(x), data.paths.door_data);
        data.Descriptors.NonSubs.pathAvgBotnecksDoor  = cellfun(@(x) mean(x), data.paths.door_data);
        % >> surface subunits
        data.Descriptors.NonSubs.edgeRGDconc       = zeros(num_2Dedgesubs, 1);
        data.Descriptors.NonSubs.edgeRGDaccessible = zeros(num_2Dedgesubs, 1);
        for i = 1 : num_2Dedgesubs
            data.Descriptors.NonSubs.edgeRGDconc(i)       = edge_2Dsubs{i}.RGDconc;
            data.Descriptors.NonSubs.edgeRGDaccessible(i) = edge_2Dsubs{i}.RGDaccessible;
        end
        % >> regions
        data.Descriptors.NonSubs.regionVols_1um     = region_vols{1};
        data.Descriptors.NonSubs.regionVols_10um    = region_vols{2};
        data.Descriptors.NonSubs.regionVols_30um    = region_vols{3};
        data.Descriptors.NonSubs.regionVols_60um    = region_vols{4};
        % Fill subunit descriptors
        data.Descriptors.Subs.volume                = zeros(num_subs, 1);
        data.Descriptors.Subs.surfArea              = zeros(num_subs, 1);
        data.Descriptors.Subs.charLength            = zeros(num_subs, 1);
        data.Descriptors.Subs.convLength            = zeros(num_subs, 1);
        data.Descriptors.Subs.avgInternalDiam       = zeros(num_subs, 1);
        data.Descriptors.Subs.aspectRatio           = zeros(num_subs, 1);
        data.Descriptors.Subs.numHalls              = zeros(num_subs, 1);
        data.Descriptors.Subs.numCrawlSpaces        = zeros(num_subs, 1);
        data.Descriptors.Subs.numConnectPores       = zeros(num_subs, 1);
        data.Descriptors.Subs.numSurroundPores      = zeros(num_subs, 1);
        data.Descriptors.Subs.numParticleEdges      = zeros(num_subs, 1);
        data.Descriptors.Subs.normNeigh             = zeros(num_subs, 1);
        data.Descriptors.Subs.meanThickness         = zeros(num_subs, 1);
        data.Descriptors.Subs.numPeaks              = zeros(num_subs, 1);
        data.Descriptors.Subs.numSurroundBeads      = zeros(num_subs, 1);
        data.Descriptors.Subs.largestSphereDiam     = zeros(num_subs, 1);
        %data.Descriptors.Subs.largestSphereVol      = zeros(num_subs, 1);
        data.Descriptors.Subs.largestDoorDiam       = zeros(num_subs, 1);
        data.Descriptors.Subs.smallestDoorDiam      = zeros(num_subs, 1);
        %data.Descriptors.Subs.smallestDoorVol       = zeros(num_subs, 1);
        data.Descriptors.Subs.axis1                 = zeros(num_subs, 1);
        data.Descriptors.Subs.axis2                 = zeros(num_subs, 1);
        data.Descriptors.Subs.axis3                 = zeros(num_subs, 1);
        data.Descriptors.Subs.ellipseClass          = zeros(num_subs, 1);
        data.Descriptors.Subs.isotropy              = zeros(num_subs, 1);
        data.Descriptors.Subs.RGDconc               = zeros(num_subs, 1);
        data.Descriptors.Subs.RGDaccessible         = zeros(num_subs, 1);
        data.Descriptors.Subs.SAparticle            = zeros(num_subs, 1);
        data.Descriptors.Subs.SAopen                = zeros(num_subs, 1);
        data.Descriptors.Subs.xCentroid             = zeros(num_subs, 1);
        data.Descriptors.Subs.yCentroid             = zeros(num_subs, 1);
        data.Descriptors.Subs.zCentroid             = zeros(num_subs, 1);
        data.Descriptors.Subs.uniqueID              = strings(num_subs, 1);
        for i = 1 : num_subs
            data.Descriptors.Subs.volume(i)         = subunits{i}.volume;
            data.Descriptors.Subs.surfArea(i)       = subunits{i}.surfArea;
            data.Descriptors.Subs.charLength(i)     = subunits{i}.charLength;
            data.Descriptors.Subs.convLength(i)     = subunits{i}.convLength;
            data.Descriptors.Subs.avgInternalDiam(i)= subunits{i}.avgInternalDiam;
            data.Descriptors.Subs.aspectRatio(i)    = subunits{i}.convLength / subunits{i}.avgInternalDiam;
            data.Descriptors.Subs.numHalls(i)       = subunits{i}.numHalls;
            data.Descriptors.Subs.numCrawlSpaces(i) = subunits{i}.numCrawlSpaces;
            data.Descriptors.Subs.numConnectPores(i)= subunits{i}.numConnectPores;
            data.Descriptors.Subs.numSurroundPores(i) = subunits{i}.numSurroundPores;
            data.Descriptors.Subs.numParticleEdges(i) = subunits{i}.numFauxEdges;
            data.Descriptors.Subs.normNeigh(i)      = subunits{i}.normNeigh;
            data.Descriptors.Subs.meanThickness(i)  = subunits{i}.meanLocalThickness;
            data.Descriptors.Subs.numPeaks(i)       = subunits{i}.numPeaks;
            data.Descriptors.Subs.numSurroundBeads(i)  = numel(subunits{i}.beadNeighbors);
            data.Descriptors.Subs.largestSphereDiam(i) = subunits{i}.largestSphereDiam;
            %data.Descriptors.Subs.largestSphereVol(i)  = subunits{i}.largestSphereVol;
            data.Descriptors.Subs.largestDoorDiam(i)   = subunits{i}.largestDoorDiam;
            data.Descriptors.Subs.smallestDoorDiam(i)  = subunits{i}.smallestDoorDiam;
            %data.Descriptors.Subs.smallestDoorVol(i)   = subunits{i}.smallestDoorVol;
            data.Descriptors.Subs.axis1(i)          = subunits{i}.ellipse_lengths(1);
            data.Descriptors.Subs.axis2(i)          = subunits{i}.ellipse_lengths(2);
            data.Descriptors.Subs.axis3(i)          = subunits{i}.ellipse_lengths(3);
            data.Descriptors.Subs.ellipseClass(i)   = subunits{i}.ellipse_class;
            data.Descriptors.Subs.isotropy(i)       = subunits{i}.isotropy;
            data.Descriptors.Subs.RGDconc(i)        = subunits{i}.RGDconc;
            data.Descriptors.Subs.RGDaccessible(i)  = subunits{i}.RGDaccessible;
            data.Descriptors.Subs.SAparticle(i)     = subunits{i}.SA_beads;
            data.Descriptors.Subs.SAopen(i)         = subunits{i}.SA_open;
            data.Descriptors.Subs.xCentroid(i)      = subunits{i}.centerCoord(1);
            data.Descriptors.Subs.yCentroid(i)      = subunits{i}.centerCoord(2);
            data.Descriptors.Subs.zCentroid(i)      = subunits{i}.centerCoord(3);
            data.Descriptors.Subs.uniqueID(i)       = subunits{i}.uniqueID;
        end

        tElapsed = toc(tStart);
        % tot_time = writeTime(tElapsed, tot_time, runtimes_file, 'Gather descriptor data:');

        time_log(timeLogIdx).Name = 'Gather descriptor data';
        time_log(timeLogIdx).Time = tElapsed;
        time_log(timeLogIdx).Units = 'sec';
        printTimeInfo(time_log(timeLogIdx));
        timeLogIdx = timeLogIdx + 1;
    end

    data.peaks           = peaks;
    data.edgePeaks       = peaks_e;
    data.ridges1D        = ridges1D;
    data.ridges2D        = ridges2D;
    data.Subunits        = subunits;
    data.Subunits_byhall = subs_byHall;
    data.subs_adjMat     = subs_adjacencysubs;
    data.edgeSubs        = edge_subs_log;
    data.voxels          = voxels;

    totalTimeElapsed = toc(totalTimeStart);
    time_log(timeLogIdx).Name = 'Total time';
    time_log(timeLogIdx).Time = totalTimeElapsed;
    time_log(timeLogIdx).Units = 'sec';
end

function [] = printTimeInfo(s)
    c = struct2cell(s);
    titleLen = 45;
    dotLen = titleLen - length(c{1});
    dots = repmat('.', [1, dotLen]);
    fprintf('%s%s %.5f %s\n', c{1}, dots, c{2}, c{3});
end
