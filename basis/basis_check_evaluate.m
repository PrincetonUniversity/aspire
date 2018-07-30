% BASIS_CHECK_EVALUATE Checks input to basis evaluate functions
%
% Usage
%    basis_check_evaluate(v, basis);
%
% Input
%    v: The coefficient vector to be evaluated in the basis.
%    basis: The basis object in which to evaluate.
%
% Description
%    Errors if the input does not correspond to the parameters of the basis
%    object.

function basis_check_evaluate(v, basis)
    if size(v, 1) ~= basis.count
        error('First dimension of v must be of size basis.count');
    end
end
