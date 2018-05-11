% PAD_SIGNAL Zero-pad a signal
%
% Usage
%    x_pad = pad_signal(x, sz);
%
% Input
%    x: The signal to be zero-padded.
%    sz: The desired size of the output signal. Note that only the first
%       numel(sz) dimensions will be zero-padded. The rest will be left
%       intact.
%
% Output
%    x_pad: The signal x, zero-padded along the first numel(sz) dimensions
%       to have size sz. The padding is done so that the original signal
%       occupies the central pixels of the padded signal. To recover the
%       original signal, use the `unpad_signal` function.

function x_pad = pad_signal(x, sz)
    [x, sz_roll] = unroll_dim(x, numel(sz)+1);

    sz0 = size(x);
    sz0 = [sz0 ones(1, numel(sz)-numel(sz0))];
    sz0 = sz0(1:numel(sz));

    if any(sz0 > sz)
        error('padding size must be larger than original signal');
    end

    x_pad = zeros([sz size(x, numel(sz)+1)]);

    sub = cell(1, numel(sz));
    for d = 1:numel(sz)
        mid(d) = ceil((sz(d)+1)/2);
        ext1(d) = ceil((sz0(d)-1)/2);
        ext2(d) = floor((sz0(d)-1)/2);

        sub{d} = mid(d)-ext1(d):mid(d)+ext2(d);
    end
    sub{numel(sz)+1} = 1:size(x, numel(sz)+1);

    sub_grid = cell(size(sub));
    [sub_grid{:}] = ndgrid(sub{:});
    ind = sub2ind(size(x_pad), sub_grid{:});

    x_pad(ind) = x;

    x_pad = roll_dim(x_pad, sz_roll);
end
