% Ellipsoid data for subunits
% The 'isotropy' value of a subunit is the dot product of the subunit's primary
% ellipsoid vector with the average vector among all subunits
 
% subunits:       complete subunit data, e.g., data.Subunits
% interior_only:  only consider interior subunits
% plot_on:        set to 1 to turn on plotting
% ellipsoid_only: set to 1 to ONLY plot ellipsoids (not subunit gridpoints
%                 as well)
% vectors_only:   set to 1 to ONLY plot eigenvectors
% ellipse_stats:  [eigenvector 1, eigenvector 2, eigenvector 3, 'isotropy']
%                 for each subunit
 
 
function ellipse_stats = ellipseData(subunits, interior_only, plot_on, ...
                                     with_subunits, vectors_only)
 
    num_subs = numel(subunits);
    ellipse_stats = zeros(num_subs, 4);
    
    ellipse_data = cell(num_subs, 2);
    if plot_on == 1
        if with_subunits == 0
            figure;
        else
            hold on
        end
        if vectors_only == 1
            color_start = 4;
            colors = distinguishable_colors(num_subs + color_start);
        end
    end
    for i = 1 : num_subs
        if interior_only
            if ~subunits{i}.edge
                [ellipse_data{i, 1}, ellipse_data{i, 2}] = ...
                          fittingEllipsoid(subunits{i}.gridpts, plot_on, with_subunits, vectors_only); % function at end of page
            end
        else
            [ellipse_data{i, 1}, ellipse_data{i, 2}] = ...
                          fittingEllipsoid(subunits{i}.gridpts, plot_on, with_subunits, vectors_only);
        end
    end
    % Get average vector
    xavg = mean(arrayfun(@(x) ellipse_data{x,1}(1,1), 1 : size(ellipse_data, 1)));
    yavg = mean(arrayfun(@(x) ellipse_data{x,1}(2,1), 1 : size(ellipse_data, 1)));
    zavg = mean(arrayfun(@(x) ellipse_data{x,1}(3,1), 1 : size(ellipse_data, 1)));
    avg_vec = [xavg, yavg, zavg] / norm([xavg, yavg, zavg]);
    % Compute 'isotropy' value, which is the dot product of the subunit's
    % primary axis against the average vector
    isotropy = zeros(num_subs, 1);
    for i = 1 : num_subs
        if vectors_only == 1
            hold on
            
            plot3([0, ellipse_data{i,1}(1,1)] * ellipse_data{i,2}(1), ...
                  [0, ellipse_data{i,1}(2,1)] * ellipse_data{i,2}(1), ...
                  [0, ellipse_data{i,1}(3,1)] * ellipse_data{i,2}(1), ...
                  '-', 'LineWidth', 2, 'Color', colors(i + color_start - 1, :));
%             plot3([0, avg_vec(1)] * ellipse_data{1,2}(1), ...
%                   [0, avg_vec(2)] * ellipse_data{1,2}(1), ...
%                   [0, avg_vec(3)] * ellipse_data{1,2}(1), ...
%                   '--', 'LineWidth', 2);
        end
        primary_unit = [ellipse_data{i,1}(1,1), ellipse_data{i,1}(2,1), ellipse_data{i,1}(3,1)];
        isotropy(i)  = dot(primary_unit, avg_vec);
        
        ellipse_stats(i, :) = [ellipse_data{i,2}(1,1), ellipse_data{i,2}(2,2), ellipse_data{i,2}(3,3), ...
                               isotropy(i)];
    end
end
 
function [V, D] = fittingEllipsoid(gridpts, plot_on, with_subunits, vectors_only)
    % Peter's help
    X = gridpts;
 
    if plot_on == 1 && with_subunits == 1
%         hold on
%         scatter3(X(:,1), X(:,2), X(:,3), '.')
%         axis equal
    end
 
    mu = mean(X);
    dataShifted = bsxfun(@minus, X, mu);
 
    % Columns of V are the eigenvectors (orientations)
    % D contains the eigenvalues (weights) in each of those directions
    [V, D] = eig(dataShifted' * dataShifted ./ size(X, 1));
    [D, perm] = sort(diag(D), 'descend');
    D = diag(D);
    V = V(:, perm);
 
    [x,y,z] = sphere(100);
    sz = size(x);
 
    % Transformation that will turn a sphere into the ellipse we want
    T = V * sqrt(D);
 
    ellip = bsxfun(@plus, T * [x(:)'; y(:)'; z(:)'], mu');
 
    x = reshape(ellip(1, :), sz(1), sz(2));
    y = reshape(ellip(2, :), sz(1), sz(2));
    z = reshape(ellip(3, :), sz(1), sz(2));
 
    if plot_on == 1 && vectors_only ~= 1
        hold on
        surf(x,y,z, 'FaceAlpha', 0.5);
        axis equal
    end
end
 
