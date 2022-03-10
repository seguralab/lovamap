% SSSR Find the Smallest Set of Smallest Rings of a graph
%
%  The calling syntax is
%
%      RESULT = SSSR(EDGES, MAXLENGTH);
%
%  where
%
%      EDGES     = N x 2 matrix of integers, where each row specifies an edge of an
%                  undirected graph. The two integers in each row are the indices
%                  of the nodes (endpoints) of the edge.
%      MAXLENGTH = maximum length of rings (loops) to be found.
%      RESULT    = M x (maxLength + 1) array of integers, where each row contains the
%                  node indices that create a ring (loop). Rows associated with loops
%                  that do not have a length of maxLength simply have trailing 0s to
%                  fill out the rest of the row. Note that the final non-zero node
%                  index in each row is the first node index, implying that if we
%                  were to traverse the nodes in order in any given row, we do
%                  indeed traverse a ring (loop).
%
%  For compatibility reasons, the while EDGES should contain integer values, the data
%  need not be of integer type. It will be converted to to an integer type internally.
%  In a smiliar fashion, even though RESULT will contain integer values, the matrix
%  will be returned as a double type.
%
%  This Matlab function invokes a MEX-file written by Ninjabyte Computing.

function [result] = sssr(edges, maxLength)
    arguments
        edges     (:, 2) int32
        maxLength (1, 1) int32 {mustBePositive}
    end

    if exist('sssr_mex', 'file')
        result = double(sssr_mex(edges, maxLength));
    else
        error('sssr:mexNotFound', 'The MEX-file sssr_mex could not be found.');
    end
end
