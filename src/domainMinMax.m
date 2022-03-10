function [a, b] = domainMinMax(beads)
    a = min(beads(:, 1:3));
    b = max(beads(:, 1:3));
end