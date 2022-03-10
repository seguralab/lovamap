% sortrows in reverse, i.e., last column to first column

function [sorted, perm] = sortrowsReverse(a)
    [d, perm] = sortrows(fliplr(a));
    sorted = fliplr(d);
end