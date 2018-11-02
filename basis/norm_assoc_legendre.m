% NORM_ASSOC_LEGENDRE Normalized associated Legendre polynomial
%
% Usage
%    npx = norm_assoc_legendre(j, m, x);
%
% Input
%    j, m: The indices of the polynomial such that j is non-negative and
%       |m| is less than or equal to j.
%    x: An array of values between -1 and +1 on which to evaluate.
%
% Output
%    npx: The value of the normalized associated Legendre polynomial of index
%       (j, m) at the points x. This function has L^2 norm 1 on the interval
%       [-1, 1]. The (unnormalized) associated Legendre `px` is related to its
%       normalized counterpart `npx` by the relation
%
%          npx = sqrt(2*factorial(j+m)/((2*j+1)*factorial(j-m)))*px;

function px = norm_assoc_legendre(j, m, x)
    % For negative m, flip sign and use the symmetry identity. In the rest, we
    % assume that m is non-negative.
    if m < 0
        m = -m;
        px = norm_assoc_legendre(j, m, x);
        px = (-1)^m*px;
        return;
    end

    % Initialize the recurrence at (m, m) and (m, m+1).
    p0 = (-1)^m*sqrt((2*m+1)/2*prod([(2*m-1):-2:1]./[2*m:-2:1]))*(1-x.^2).^(m/2);

    p1 = x.*(sqrt(2*m+3)*p0);

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
        px = sqrt((2*n+3)/((n+1+m)*(n+1-m)))* ...
            (sqrt(2*n+1)*x.*p1-sqrt((n+m)*(n-m)/(2*n-1))*p0);

        p0 = p1;
        p1 = px;
    end
end
