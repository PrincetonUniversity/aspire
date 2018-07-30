% VEC_TO_IM Unroll vectors to images
%
% Usage
%    im = vec_to_im(vec);
%
% Input
%    vec: An N^2-by-... array.
%
% Output
%    im: An N-by-N-by-... array.

function im = vec_to_im(vec)
    N = round(size(vec, 1)^(1/2));

    if size(vec, 1) ~= N^2
        error('images must be square');
    end

    sz = size(vec);
    im = reshape(vec, [N*ones(1, 2) sz(2:end)]);
end
