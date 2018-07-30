% ROLL_DIM Roll trailing dimensions back up
%
% Usage
%    X = roll_dim(X, roll_sz);
%
% Input
%    X: The array whose trailing dimensions are to be rolled back up.
%    roll_sz: The original size of the trailing dimensions of X before being
%       unrolled.
%
% Output
%    X: The original signal X with the last dimension reshaped to have size
%       roll_sz.

function X = roll_dim(X, roll_sz)
    if length(roll_sz) > 1
        sz = size(X);
        if sz(2) == 1
            sz = sz(1);
        end

        X = reshape(X, [sz(1:end-1) roll_sz]);
    end
end
