% SUBMATFUN Apply function to submatrices
%
% Usage
%    out = submatfun(fun, in, dim);
%
% Input
%    fun: The function handle to apply to each submatrix.
%    in: The input matrix.
%    dim: The dimension along which to index the submatrices.
%
% Output
%    out: The matrix consisting of the processed submatrices.

function out = submatfun(fun, in, dim)
    % Get original dimensions
    in_sz = size(in);
    in_dims = ndims(in);

    if dim > in_dims
        in_sz = [in_sz ones(1, dim-in_dims)];
        in_dims = dim;
    end

    % Extract first slice and process
    cln = {':'};
    idx.type = '()';
    idx.subs = [cln(ones(1,dim-1)), 1, cln(ones(1,in_dims-dim))];
    in_slice = permute(subsref(in, idx), [1:dim-1 dim+1:in_dims dim]);
    out_slice = fun(in_slice);
    out_slice = permute(out_slice, [1:dim-1 in_dims dim:in_dims-1]);

    % Determine size of output
    out_sz = size(out_slice);
    out_sz = [out_sz ones(1, in_dims-numel(out_sz))];
    out_sz(dim) = in_sz(dim);

    % Create output array
    out = zeros(out_sz, class(out_slice));

    % Assign first slice output
    out = subsasgn(out, idx, out_slice);

    % Process remaining slices
    for k = 2:in_sz(dim)
        idx.subs = [cln(ones(1,dim-1)), k, cln(ones(1,in_dims-dim))];
        in_slice = permute(subsref(in, idx), [1:dim-1 dim+1:in_dims dim]);
        out_slice = fun(in_slice);
        out_slice = permute(out_slice, [1:dim-1 in_dims dim:in_dims-1]);
        out = subsasgn(out, idx, out_slice);
    end
end
