% MAT_TO_BLK_DIAG Convert full matrix to block diagonal matrix
%
% Usage
%    blk_diag = mat_to_blk_diag(mat, blk_partition);
%
% Input
%    mat: A full matrix.
%    blk_partition: A matrix block partition in the form of a K-by-2 array,
%       where `blk_partition(:,1)` corresponds to a partition of the rows and
%       `blk_partition(:,2)` is a partition of the columns.
%
% Output
%    blk_diag: A block diagonal matrix consisting of the `K` diagonal blocks
%       of `mat`. This is in the form of a cell array containing `K` elements,
%       each a separate block.
%
% See also
%    blk_diag_to_mat

function blk_diag = mat_to_blk_diag(mat, blk_partition)
    rows = blk_partition(:,1);
    cols = blk_partition(:,2);

    blk = mat2cell(mat, rows, cols);

    blk_diag = diag(blk);
end
