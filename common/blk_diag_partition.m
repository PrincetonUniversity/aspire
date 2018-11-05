% BLK_DIAG_PARTITION Get partition of block diagonal matrix
%
% Usage
%    blk_partition = blk_diag_partition(blk_diag);
%
% Input
%    blk_diag: A block diagonal matrix in the form of a cell array. Each
%       element corresponds to a diagonal block.
%
% Output
%    blk_partition: The matrix block partition of `blk_diag` in the form of a
%       K-by-2 array, where `blk_partition(:,1)` corresponds to a partition of
%       the rows and `blk_partition(:,2)` is a partition of the columns.

function blk_partition = blk_diag_partition(blk_diag)
    blk_diag = blk_diag(:);

    blk_partition = cellfun(@size, blk_diag, 'uniformoutput', false);

    blk_partition = cat(1, blk_partition{:});
end
