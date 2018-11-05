% BLK_DIAG_ZEROS Construct a block diagonal zero matrix
%
% Usage
%    blk_diag = blk_diag_zeros(blk_partition, precision);
%
% Input
%    blk_partition: A matrix block partition in the form of a K-by-2 array,
%       where `blk_partition(:,1)` corresponds to a partition of the rows and
%       `blk_partition(:,2)` is a partition of the columns.
%    precision: The precision, 'single' or 'double', of the matrix (default
%       'double').
%
% Output
%    blk_diag: A block diagonal matrix consisting of `K` zero blocks.

function blk_diag = blk_diag_zeros(blk_partition, precision)
    if nargin < 2 || isempty(precision)
        precision = 'double';
    end

    blk_partition = num2cell(blk_partition, 2);

    blk_diag = cellfun(@(sz)(zeros(sz, precision)), blk_partition, ...
        'uniformoutput', false);

    blk_diag = blk_diag(:);
end
