% POWER_METHOD    Iterative method for computing the largest eigenvalue, in
%                 magnitude, and its coresponding eigenvector of a
%                 real symmetric matrix.
%
%   [LAMBDA, V] = POWER_METHOD(A, NUM_ITERS, TOL) returns the largest
%       eigenvalue, LAMBDA, in magnitude, and its corresponding
%       eigenvector, V. This functions requires the following as arguments:
%
%       A         - Real symmetric matrix (necessarily square)
%       NUM_ITERS - Maximum number of iterations for the power method.
%       TOL       - The power method will exit early if the absolute value
%                   of the difference between two successive approximations
%                   of the eigenvalue is less than TOL.
%
%    Important: This method randomly initializes its starting eigenvector
%    approximation (for good reason), which implies subsequent runs of this 
%    method using the same input arguments will likely yield different 
%    results. However, if the number of iterations is sufficiently high and 
%    the tolerance is sufficiently small, then subsequent results should 
%    only differ by a trivial amount.
%
% Author: Petes

function [lambda_k1, vk] = power_method(A, numIters, tol)

    % Matrix must be square
    [nr, nc] = size(A);
    if nr ~= nc
        error('Matrix must be square.');
    end

    % Randomly generate an initial vector and normalize it.
    % We randomly generate this vector so, in general, it has a higher
    % probability not being perpendicular to the eigenvector of largest
    % magnitude.
    vk = rand(nr, 1);
    vk = vk / norm(vk);

    % Perform one iteration first to initialize an approximation for the
    % eignevalue. Note that, while we use the Rayleigh quotient to estimate the
    % eignevalue, we do not need to divide by vk' * vk since we're always
    % working with normalized vectors.
    vk1 = A * vk;
    lambda_k1 = dot(vk, vk1);
    vk = vk1 / norm(vk1);

    % Iterate until max number of iterations or if the difference between two
    % consecutive eigenvalue approximations is small enough.
    for it = 2 : numIters
        lambda_k = lambda_k1;

        vk1 = A * vk;
        lambda_k1 = dot(vk, vk1);
        vk = vk1 / norm(vk1);

        if abs(lambda_k1 - lambda_k) < tol
            return
        end 
    end

end
