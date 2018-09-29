% BLK_DIAG_ADD Add block diagonal matrices
%
% Usage
%    blk_diag_C = blk_diag_add(blk_diag_A, blk_diag_B);
%
% Input
%    blk_diag_A, blk_diag_B: Two block diagonal matrices in the form of cell
%       arrays. Each element corresponds to a diagonal block.
%
% Output
%    blk_diag_C: The block diagonal array corresponding to the sum of the two
%       input matrices.

function blk_diag_C = blk_diag_add(blk_diag_A, blk_diag_B)
    blk_diag_A = blk_diag_A(:);
    blk_diag_B = blk_diag_B(:);

    compat_fun = @(block1, block2)(all(size(block1) == size(block2)));
    if ~all(cellfun(compat_fun, blk_diag_A, blk_diag_B))
        error('Block diagonal matrix sizes are incompatible.');
    end

    blk_diag_C = cellfun(@plus, blk_diag_A, blk_diag_B, ...
        'uniformoutput', false);
end
