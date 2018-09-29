% BLK_DIAG_APPLY Multiply a block diagonal matrix by a vector
%
% Usage
%    y = blk_diag_apply(blk_diag, x);
%
% Input
%    blk_diag: A block diagonal matrix in the form of a cell array. Each
%       element corresponds to a diagonal block.
%    x: A vector to be multiplied by `blk_diag`. May also be a matrix, in
%       which case each column is multiplied separately.
%
% Output
%    y: The result of left-multiplying `x` by `blk_diag`. If `x` is a matrix,
%       `y` is also a matrix.
%
% See also
%    blk_diag_solve

function y = blk_diag_apply(blk_diag, x)
    blk_diag = blk_diag(:);

    cols = cellfun(@(block)(size(block, 2)), blk_diag);

    if sum(cols) ~= size(x, 1)
        error('Sizes of matrix `blk_diag` and `x` are not compatible.');
    end

    x = mat2cell(x, cols, size(x, 2));

    y = cellfun(@mtimes, blk_diag, x, 'uniformoutput', false);

    y = cat(1, y{:});
end
