% Label 3D axes
%
% varargin : set to 1 if you want to shift the graph to the origin


function axes_labels(fig, fig_title, font_size, varargin)

    title(fig_title, 'FontSize', font_size);
    xlabel('x (\mu\itm)', 'FontSize', font_size);
    ylabel('y (\mu\itm)', 'FontSize', font_size);
    zlabel('z (\mu\itm)', 'FontSize', font_size);

    xt    = get(gca(fig), 'XTick');
    xtlab = get(gca(fig), 'XTickLabel');
    yt    = get(gca(fig), 'YTick');
    ytlab = get(gca(fig), 'YTickLabel');
    if ~isempty(varargin)
        if xt(1) == 0
            set(gca(fig), 'YTick', xt, 'YTickLabel', xtlab, 'FontSize', 20);
            set(gca(fig), 'XTick', xt, 'XTickLabel', xtlab, 'FontSize', 20);
        elseif yt(1) == 0
            set(gca(fig), 'XTick', yt, 'XTickLabel', ytlab, 'FontSize', 20);
            set(gca(fig), 'YTick', yt, 'YTickLabel', ytlab, 'FontSize', 20);
        else
            tickspace_x = 0; %tickspace_x = floor(abs(xt(2) - xt(1)));
            set(gca(fig), 'XTick', xt, 'XTickLabel', ...
                (xt - xt(1)) + tickspace_x, 'FontSize', 20);
            set(gca(fig), 'YTick', yt, 'YTickLabel', ...
                (yt - yt(1)) + tickspace_x, 'FontSize', 20);
        end
        zt    = get(gca(fig), 'ZTick');
        tickspace_z = 0; %tickspace_z = floor(abs(zt(2) - zt(1)));
        if zt(1) ~= 0
            set(gca(fig), 'ZTick', zt, 'ZTickLabel', ...
            (zt - zt(1)) + tickspace_z, 'FontSize', 20);
        else
            set(gca(fig), 'FontSize', 20);
        end
    else
        set(gca(fig), 'FontSize', 20);
    end
end