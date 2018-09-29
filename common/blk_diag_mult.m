% BLK_DIAG_MULT Multiply block diagonal matrices
%
% Usage
%    blk_diag_C = blk_diag_add(blk_diag_A, blk_diag_B);
%
% Input
%    blk_diag_A, blk_diag_B: Two block diagonal matrices in the form of cell
%       arrays. Each element corresponds to a diagonal block.
%
% Output
%    blk_diag_C: The block diagonal array corresponding to the product of the
%       two input matrices.

function blk_diag_C = blk_diag_mult(blk_diag_A, blk_diag_B)
    blk_diag_A = blk_diag_A(:);
    blk_diag_B = blk_diag_B(:);

    if isnumeric(blk_diag_B)
        blk_diag_C = blk_diag_mult(blk_diag_B, blk_diag_A);
        return;
    end

    if isnumeric(blk_diag_A)
        blk_diag_C = cellfun(@(block)(blk_diag_A*block), blk_diag_B, ...
            'uniformoutput', false);
        return;
    end

    compat_fun = @(block1, block2)(size(block1, 2) == size(block2, 1));
    if ~all(cellfun(compat_fun, blk_diag_A, blk_diag_B))
        error('Block diagonal matrix sizes are incompatible.');
    end

    blk_diag_C = cellfun(@mtimes, blk_diag_A, blk_diag_B, ...
        'uniformoutput', false);
end
