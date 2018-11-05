% BLK_DIAG_TO_MAT Convert block diagonal matrix to full matrix
%
% Usage
%    mat = blk_diag_to_mat(blk_diag);
%
% Input
%    blk_diag: A block diagonal matrix in the form of a cell array. Each
%       element corresponds to a diagonal block.
%
% Output
%    mat: The full matrix corresponding to the block diagonal matrix
%       `blk_diag`.
%
% See also
%    mat_to_blk_diag

function mat = blk_diag_to_mat(blk_diag)
    mat = blkdiag(blk_diag{:});
end
