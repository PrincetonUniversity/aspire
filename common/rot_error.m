% ROT_ERR Calculates the MSE between two arrays of rotation matrices
%
% Usage
%    [mse, Q] = rot_error(rots1, rots2);
%
% Input
%    rots1, rots2: Two 3-by-3-by-n arrays of rotation matrices, indexed along
%       the third dimension.
%
% Output
%    mse: The minimum of
%
%          1/n*sum_{k=1}^n |Q*rots1(:,:,k)-rots2(:,:,k)|_F^2 ,
%
%       where Q ranges over the set of orthogonal matrices.
%    Q: The 3-by-3 orthogonal matrix Q that minimizes the above objective.

function [mse, Q] = rot_error(rots1, rots2)
    if size(rots1, 3) ~= size(rots2, 3)
        error('Rotation arrays must have the same number of matrices.');
    end

    A = zeros(3);

    for k = 1:size(rots1, 3)
        A = A + rots2(:,:,k)*rots1(:,:,k)';
    end

    Q = closest_orthogonal_matrix(A);

    rots1 = submatfun(@(X)(Q*X), rots1, 3);

    mse = tnorm(rots1-rots2)^2/size(rots1, 3);
end
