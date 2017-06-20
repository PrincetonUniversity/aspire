% UNPAD_SIGNAL Extract central part of signal
%
% Usage
%    x = unpad_signal(x_pad, sz);
%
% Input
%    x_pad: The signal whose center should be extracted.
%    sz: The size of the central box to extract. Note that only the first
%       numel(sz) dimensions of x will be indexed in this extraction. The
%       remaining dimensions will be untouched.
%
% Output
%    x: The central box of size sz extracted from the first numel(sz)
%       dimensions of x. To recover a signal of size x with the discarded
%       parts zeroed out, use the `pad_signal` function.

function x = unpad_signal(x_pad, sz)
    [x_pad, sz_roll] = unroll_dim(x_pad, numel(sz)+1);

    % TODO: The similarities with `pad_signal` are significant. Needs to be
    % reduced.
    sz0 = size(x_pad);
    sz0 = [sz0 ones(1, numel(sz)-numel(sz0))];
    sz0 = sz0(1:numel(sz));

    if any(sz0 < sz)
        error('extraction box must be smaller than original signal')
    end

    x = zeros([sz size(x_pad, numel(sz)+1)]);

    sub = cell(1, numel(sz)+1);
    for d = 1:numel(sz)
        mid(d) = ceil((sz0(d)+1)/2);
        ext1(d) = ceil((sz(d)-1)/2);
        ext2(d) = floor((sz(d)-1)/2);

        sub{d} = mid(d)-ext1(d):mid(d)+ext2(d);
    end
    sub{numel(sz)+1} = 1:size(x_pad, numel(sz)+1);

    sub_grid = cell(size(sub));
    [sub_grid{:}] = ndgrid(sub{:});
    ind = sub2ind(size(x_pad), sub_grid{:});

    x = x_pad(ind);

    x = roll_dim(x, sz_roll);
end
