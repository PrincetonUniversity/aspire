% ROLL_DIM Roll trailing dimensions back up
%
% Usage
%    X = roll_dim(X, sz_roll);
%
% Input
%    X: The array whose trailing dimensions are to be rolled back up.
%    sz_roll: The original size of the trailing dimensions of X before being
%       unrolled.
%
% Output
%    X: The original signal X with the last dimension reshaped to have size
%       sz_roll.

function X = roll_dim(X, sz_roll)
    if length(sz_roll) > 1
        sz = size(X);
        if sz(2) == 1
            sz = sz(1);
        end

        X = reshape(X, [sz(1:end-1) sz_roll]);
    end
end
