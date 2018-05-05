% IM_TO_VEC Roll up images into vectors
%
% Usage
%    vec = im_to_vec(im);
%
% Input
%    im: An N-by-N-by-... array.
%
% Output
%    vec: An N^2-by-... array.

function vec = im_to_vec(im)
    N = size(im, 1);

    if size(im, 2) ~= N
        error('Images in `im` must be square.');
    end

    sz = size(im);
    vec = reshape(im, [N^2 sz(3:end) 1]);
end
