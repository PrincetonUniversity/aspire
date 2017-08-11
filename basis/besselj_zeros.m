% BESSELJ_ZEROS Find k first zeros of J_nu
%
% Usage
%    z = besselj_zeros(nu, k);
%
% Input
%    nu: Order of the Bessel function J_nu. Must be less than 1e7.
%    k: The number of zeros desired.
%
% Output
%    z: The first k zeros of J_nu, ordered by increasing magnitude.

% Adapted from "zerobess.m" by Jonas Lundgren <splinefit@gmail.com>.

function z = besselj_zeros(nu, k)
    if nu < 0 || nu > 1e7
        error('Input `nu` must be between 0 and 1e7.');
    end

    z = zeros(k, 1);

    % Guess first zeros using powers of nu.
    c0 = [ 0.1701 -0.6563 1.0355 1.8558; ...
           0.1608 -1.0189 3.1348 3.2447; ...
          -0.2005 -1.2542 5.7249 4.3817];

    z0 = nu + c0*((nu+1).^[-1; -2/3; -1/3; 1/3]);

    % Refine guesses.
    z(1:3) = besselj_newton(nu, z0);

    n = 3;
    j = 2;
    err_tol = 5e-3;

    % Estimate further zeros iteratively using spacing of last three zeros so
    % far.
    while n < k
        j = min(j, k-n);

        % Use last three zeros to predict spacing for next j zeros.
        r = diff(z(n-2:n)) - pi;
        if r(1)*r(2) > 0 && r(1)/r(2) > 1
            p = log(r(1)/r(2))/log(1-1/(n-1));
            dz = pi + r(2)*exp(p*log(1+[1:j]'/(n-1)));
        else
            dz = pi*ones(j, 1);
        end

        % Guess and refine.
        z0 = z(n) + cumsum(dz);
        z(n+1:n+j) = besselj_newton(nu, z0);

        % Check to see that the sequence of zeros makes sense.
        if ~check_besselj_zeros(nu, z(n-1:n+j))
            error('Unable to properly estimate Bessel function zeros.');
        end

        % Check how far off we are.
        err = (z(n+1:n+j) - z0)./diff(z(n:n+j));

        n = n+j;

        if max(abs(err)) < err_tol
            % Predictions were close enough; double number of zeros.
            j = 2*j;
        else
            % Some prediction were off; set to double the number of good
            % predictions.
            j = 2*find(abs(err) >= err_tol, 1);
        end
    end

    z = z';
end

function b = check_besselj_zeros(nu, z)
    dz = diff(z);
    ddz = diff(dz);

    b = true;

    b = b && all(isreal(z));
    b = b && z(1) > 0;
    b = b && all(dz > 3);

    if nu >= 0.5
        b = b && all(ddz < 16*eps(z(2:end-1)));
    else
        b = b && all(ddz > -16*eps(z(2:end-1)));
    end
end
