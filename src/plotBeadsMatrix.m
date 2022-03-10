% Plot beads from the labeled domain matrix
% beads:            cell array of bead data
% voxels:           [x, y, z] grid data
% bead_diams:       data.Global.beadDiams
% bead_color:       'r', 'g', 'b', or 'y' - to set all beads to this color
%                   'auto' - for auto coloring based on beadColor.m
% varargin:         either: 1) cell array of subunits implying only beads touching these
%                   subunits should be plotted, or 2) 'all' indicating plot all
%                   beads

function [] = plotBeadsMatrix(f, beads, bead_diams, bead_color, voxels, ...
                              domain, bead_alpha, dateStamp, gif, varargin)

    figure(f);
    hold on
    
    % Grab beads to plot
    if ~isempty(varargin) && iscell(varargin{1})
        % want neighboring beads of specific subunits
        subunits = varargin{1};
        bead_cache = [];
        bead_cachelog = false(length(beads), 1);
        for i = 1 : length(subunits)
            bead_cache = [bead_cache; subunits{i}.beadNeighbors];
        end
        bead_cache = unique(bead_cache);
        bead_cachelog(bead_cache) = true;
    elseif ~isempty(varargin) && strcmp(varargin{1}, 'all')
        % plot all beads
        bead_cachelog = true(length(beads), 1);
    else
        % plot all beads
        bead_cachelog = true(length(beads), 1);
        xlim([domain(1) domain(2)]);
        ylim([domain(3) domain(4)]);
        zlim([domain(5) domain(6)]);
        axes_labels(f, '', 20);
    end
    beads2plot = beads(bead_cachelog);
    num_beads = length(beads2plot);
    
    switch bead_color
        case 'auto'
            % Color by particle size
            % round volumes to nearest 10
            diams = round(bead_diams / 10) * 10;
        case 'r'
            bcolor = '#C60806';
        case 'g'
            bcolor = [ 76, 187,  23] / 255;
        case 'b'
            bcolor = [ 15,  82, 186] / 255;
        case 'y'
            bcolor = [255, 242,   0] / 255;
        otherwise
            bcolor = bead_color;
    end

    for i = 1 : num_beads
        x = voxels(beads2plot{i}, 1);
        y = voxels(beads2plot{i}, 2);
        z = voxels(beads2plot{i}, 3);
        k = boundary(x, y, z);

        hold on

        if strcmp(bead_color, 'auto')
            trisurf(k, x, y, z, 'Facecolor', beadColor(diams(i)), 'FaceAlpha', 1, ...
                    'EdgeColor', 'none');
        else
            trisurf(k, x, y, z, 'Facecolor', bcolor, 'FaceAlpha', 1, ...
                    'EdgeColor', 'none');
        end
    end
    alpha(bead_alpha);

    if gif == 1
        fig_rot_animation(f, ['./gifs/ParticleDomain_', num2str(dateStamp), '.gif']);
    end
end