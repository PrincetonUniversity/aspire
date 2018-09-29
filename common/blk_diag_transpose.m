% BLK_DIAG_TRANSPOSE Transpose a block diagonal matrix
%
% Usage
%    blk_diag = blk_diag_transpose(blk_diag);
%
% Input
%    blk_diag: A block diagonal matrix in the form of a cell array. Each
%       element corresponds to a diagonal block.
%
% Output
%    blk_diag: The same block diagonal matrix with each block transposed.

function blk_diag = blk_diag_transpose(blk_diag)
    blk_diag = blk_diag(:);

    blk_diag = cellfun(@transpose, blk_diag, 'uniformoutput', false);
end
