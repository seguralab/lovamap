% Plot 3D beads

function [] = plotBeads3D(bead_data, bead_alpha, sphere_res)

    % Generate (x,y,z) coordinates of a unit with sphere_res resolution
    [x, y, z] = sphere(sphere_res);
    
    % Plot beads
    for i = 1 : size(bead_data, 1)
        if bead_data(i, 4) > 1000 % adjust this when creating color threshold
            surfl(x * bead_data(i,4) + bead_data(i,1), ...
                  y * bead_data(i,4) + bead_data(i,2), ...
                  z * bead_data(i,4) + bead_data(i,3));
            hold on
            material dull
            shading interp      

            hold on
            
            % Blue particles
            vec_length = 100;
            dkblue = [24, 60, 148]/255;
            blue = [15, 82, 186]/255;
            color_blue = [linspace(dkblue(1), blue(1), vec_length)', ...
                          linspace(dkblue(2), blue(2), vec_length)', ...
                          linspace(dkblue(3), blue(3), vec_length)'];
            colormap(color_blue)
        else
            surfl(x * bead_data(i,4) + bead_data(i,1), ...
                  y * bead_data(i,4) + bead_data(i,2), ...
                  z * bead_data(i,4) + bead_data(i,3));
            hold on
            material dull
            shading interp      

            hold on
            
            % Red particles
            vec_length = 100;
            dkred = [95, 2, 8]/255;
            red = [198, 8, 6]/255;
            color_red = [linspace(dkred(1), red(1), vec_length)', ...
                         linspace(dkred(2), red(2), vec_length)', ...
                         linspace(dkred(3), red(3), vec_length)'];
            colormap(color_red)
        end
    end

    grid on
    alpha(bead_alpha)
    hold off
    axis equal
end

% set(gca, 'Color', 'k');
%     color = get(f, 'Color');
%     set(gca, 'XColor', color, 'YColor', color, 'TickDir', 'out');
    
% cmap = colormap(f);
% int_vs = integrin_map(integrin_map >= 0);
% cbins = linspace(min(int_vs(:)), max(integrin_map(:)), size(cmap, 1) + 1);
% alpha_vals = linspace(0, 1, size(cmap, 1));



% % For plotting a cylinder
% Radius_Container = 1;
% [XX1, YY1, ZZ1] = cylinder(Radius_Container,200);
% surf(XX1 + Radius_Container, YY1 + Radius_Container, ZZ1)