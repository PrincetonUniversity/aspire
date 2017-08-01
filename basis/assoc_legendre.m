% ASSOC_LEGENDRE Associated Legendre polynomial
%
% Usage
%    px = assoc_legendre(j, m, x);
%
% Input
%    j, m: The indices of the polynomial such that j is non-negative and
%       |m| is less than or equal to j.
%    x: An array of values between -1 and +1 on which to evaluate.
%
% Output
%    px: The value of the associated Legendre polynomial of index (j, m) at
%       the points x.

function px = assoc_legendre(j, m, x)
    % For negative m, flip sign and use the symmetry identity. In the rest, we
    % assume that m is non-negative.
    if m < 0
        m = -m;
        px = assoc_legendre(j, m, x);
        px = (-1)^m*1/prod((j-m+1):(j+m))*px;
        return;
    end

    % Initialize the recurrence at (m, m) and (m, m+1).
    p0 = (-1)^m*prod((2*m-1):-2:1)*(1-x.^2).^(m/2);
    p1 = x.*((2*m+1)*p0);

    % If these are the desired indices, return these initial values.
    if j == m
        px = p0;
        return;
    elseif j == m+1
        px = p1;
        return;
    end

    % Fixing m, work our way up from (m, m+1) to (m, j).
    for n = m+1:j-1
        px = 1/(n-m+1)*((2*n+1)*x.*p1-(n+m)*p0);

        p0 = p1;
        p1 = px;
    end
end
