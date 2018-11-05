% RADIAL_FUN_TO_IMAGE Evaluate a radial function over a 2D grid
%
% Usage
%    im = radial_fun_to_image(fun, L);
%
% Input
%    fun: A function handle taking one input: a radius between 0
%       and sqrt(2).
%
% Output
%    im: The function values on a square grid spanning [-1/2, 1/2] by [-1/2,
%        1/2].

function im = radial_fun_to_image(fun, L)
    m2d = mesh_2d(L);

    im = fun(m2d.r/2);
end
