function [out] = all_sub2ind(res, idx)
% ALL_SUB2IND Helper function for calculating the linear indices for a
% list of 3D integer indices given a fixed resolution.
%
%   OUT = ALL_SUB2IND(RES, IDX) where
%   RES = size of the matrix, must be 1 x 3
%   IDX = N x 3 array of all 3D subscripts that need to be
%         converted to their associated linear indices

    % Peter is amazeballs

    % Early out if there are no subscripts to convert
    if isempty(idx)
        out = [];
    end

    res = cast(res(:)', class(idx));

    dim = length(res);
    strides = repmat(cumprod(res(1 : end - 1)), size(idx, 1), 1);
    out = idx(:, 1);
    for j = 2 : dim
        out = out + (idx(:, j) - 1) .* strides(:, j - 1);
    end
end
