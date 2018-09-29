% BLK_DIAG_EYE Construct a block diagonal identity matrix
%
% Usage
%    blk_diag = blk_diag_eye(blk_partition, precision);
%
% Input
%    blk_partition: A matrix block partition in the form of a K-by-2 array,
%       where `blk_partition(:,1)` corresponds to a partition of the rows and
%       `blk_partition(:,2)` is a partition of the columns.
%    precision: The precision, 'single' or 'double', of the matrix (default
%       'double').
%
% Output
%    blk_diag: A block diagonal matrix consisting of `K` identity blocks.

function blk_diag = blk_diag_eye(blk_partition, precision)
    if nargin < 2 || isempty(precision)
        precision = 'double';
    end

    blk_partition = num2cell(blk_partition, 2);

    blk_diag = cellfun(@(sz)(eye(sz, precision)), blk_partition, ...
        'uniformoutput', false);

    blk_diag = blk_diag(:);
end
