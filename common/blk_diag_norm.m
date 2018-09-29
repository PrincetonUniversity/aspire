% BLK_DIAG_NORM Compute the 2-norm of a block diagonal matrix
%
% Usage
%    n = blk_diag_norm(blk_diag);
%
% Input
%    blk_diag: A block diagonal matrix in the form of a cell array. Each array
%       corresponds to a diagonal block.
%
% Output
%    n: The 2-norm of the block diagonal matrix.

function n = blk_diag_norm(blk_diag)
    n = max(cellfun(@norm, blk_diag));
end
