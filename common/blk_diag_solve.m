% BLK_DIAG_SOLVE Solve a linear system involving a block diagonal matrix
%
% Usage
%    x = blk_diag_solve(blk_diag, y);
%
% Input
%    blk_diag: A block diagonal matrix in the form of a cell array. Each
%       element corresponds to a diagonal block.
%    y: The right-hand side in the linear system.  May also be a matrix, in
%       which case each column is solved for separately.
%
% Output
%    y: The result of solving the linear system formed by the matrix
%       `blk_diag` and the right-hand size `y`. If `y` is a matrix, `x`
%       is also a matrix.
%
% See also
%    blk_diag_apply

function x = blk_diag_solve(blk_diag, y)
    blk_diag = blk_diag(:);

    rows = cellfun(@(block)(size(block, 1)), blk_diag);

    if sum(rows) ~= size(y, 1)
        error('Sizes of matrix `blk_diag` and `x` are not compatible.');
    end

    y = mat2cell(y, rows, size(y, 2));

    x = cellfun(@mldivide, blk_diag, y, 'uniformoutput', false);

    x = cat(1, x{:});
end
