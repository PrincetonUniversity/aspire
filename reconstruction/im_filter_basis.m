% IM_FILTER_BASIS Matrix representation of image filtering in basis
%
% Usage
%    H = im_filter_basis(filter_f, basis);
%
% Input
%    filter_f: An array of filters in the Fourier domain (centered) of size
%       basis.sz-by-K, where K is the number of filters.
%    basis: A basis object.
%
% Output
%    H: An array of matrices representing the application of the filters in the
%       basis. It is of size basis.count-by-basis.count-by-K.

function H = im_filter_basis(filter_f, basis)
    filter_ct = size(filter_f, 3);

    H = zeros([basis.count*ones(1, 2) filter_ct], class(filter_f));

    for k = 1:filter_ct
        fun = @(im)(real(im_filter(im, filter_f(:,:,k))));

        H(:,:,k) = proj_fun_basis(fun, basis);
    end
end
