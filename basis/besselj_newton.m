% BESSELJ_NEWTON Find closest zero using Newton's method
%
% Usage
%    z = besselj_newton(nu, z0, max_iter);
%
% Input
%    nu: The order of the Bessel function.
%    z0: The initial value at which to start the search.
%    max_iter: Maximum number of iterations (default 10).
%
% Output
%    z: The zero of J_nu closest to z0 after max_iter iterations.

function z = besselj_newton(nu, z0, max_iter)
    if nargin < 3 || isempty(max_iter)
        max_iter = 10;
    end

    z = z0;

    % Factor worse than machine precision.
    c = 8;

    for iter = 1:max_iter
        % Calculate values and derivatives at z.
        f = besselj(nu, z);
        fp = besselj(nu-1, z) - nu*f./z;

        % Update zeroes.
        dz = -f./fp;
        z = z + dz;

        % Check for convergence.
        if all(abs(dz) < c*eps(z))
            break;
        end

        % If we're not converging yet, start relaxing convergence criterion.
        if iter >= 6
            c = 2*c;
        end
    end
end
