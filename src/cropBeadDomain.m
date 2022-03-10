% Crop packed spheres

function [beads_cropped, domain_cropped] = cropBeadDomain(bead_data, domain, percentDomain)

    % Determine domain to plot
    if (domain(6) - domain(5)) <= ...
       ((domain(2) - domain(1)) * percentDomain)

        portionTrim = floor([floor((domain(2) - domain(1)) * ...
                                                    (1 - percentDomain)), ...
                       floor((domain(4) - domain(3)) * ...
                                                    (1 - percentDomain))] ...
                       ./ 2);

        xmin = domain(1) + portionTrim(1);
        xmax = domain(2) - portionTrim(1);
        ymin = domain(3) + portionTrim(2);
        ymax = domain(4) - portionTrim(2);
        zmin = domain(5);
        zmax = domain(6);
        domain_cropped = [xmin, xmax, ymin, ymax, zmin, zmax];

    else
        portionTrim = floor([floor((domain(2) - domain(1)) * ...
                                                    (1 - percentDomain)), ...
                       floor((domain(4) - domain(3)) * ...
                                                    (1 - percentDomain)), ...
                       floor((domain(6) - domain(5)) * ...
                                                    (1 - percentDomain))] ...
                       ./ 2);

        xmin = domain(1) + portionTrim(1);
        xmax = domain(2) - portionTrim(1);
        ymin = domain(3) + portionTrim(2);
        ymax = domain(4) - portionTrim(2);
        zmin = domain(5) + portionTrim(3);
        zmax = domain(6) - portionTrim(3);
        domain_cropped = [xmin, xmax, ymin, ymax, zmin, zmax];
    end

    % Store beads that lie within domain
    beads_cropped = [];
    for i = 1 : size(bead_data, 1)
        xyz = bead_data(i, 1:3);
        if (all(xyz >= [xmin, ymin, zmin]) & all(xyz <= [xmax, ymax, zmax]))
            beads_cropped = vertcat(beads_cropped, bead_data(i, :));
        end
    end
    % Store boundary data
    if isempty(beads_cropped)
        error('Must increase crop percent to include at least 1 bead')
    end
end
