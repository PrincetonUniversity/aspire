% FUN_ONTO_BASIS Project linear mapping onto basis
%
% Usage
%    A = fun_onto_basis(fun, basis_from, basis_to);
%
% Input
%    fun: A function handle mapping multidimensional vectors from one space
%       linearly to another. It must take as input arrays of the form
%       basis_from.sz-by-K and output arrays of the form basis_to.sz-by-K.
%    basis_from: The basis of the domain space to project onto.
%    basis_to: The basis of the target space to project onto (default same as 
%       basis_from).
%
% Output
%    A: A matrix of size basis_to.count-by-basis_from.count representing the
%       linear mapping fun in the bases basis_from and basis_to. In other words,
%       given a vector x of length basis_from.count, we have
%
%          A*x = basis_to.expand(fun(basis_from.evaluate(x)));

function A = fun_onto_basis(fun, basis_from, basis_to)
    if nargin < 3 || isempty(basis_to)
        basis_to = basis_from;
    end

    X = basis_from.evaluate(eye(basis_from.count));

    Y = fun(X);

    A = basis_to.expand(Y);
end
