% RADIAL_FILTER_TO_FB_MAT Represent radial filtering as Fourier-Bessel matrix
%
% Usage
%    h_fb = radial_filter_to_fb_mat(h_fun, basis);
%
% Input
%    h_fun: A function handle taking one input: a radial frequency between 0
%       and 1/2, where 1/2 is the Nyquist frequency of the signal.
%    basis: The basis in which the coefficients are expanded. Must be a
%       Fourier-Bessel basis obtained from the `fb_basis` function.
%
% Output
%    h_fb: A block diagonal matrix (may be manipulated with by the
%       `blk_diag_*` functions) representing the application of the transfer
%       function `h_fun` to an image represented in a Fourier-Bessel basis.
%
% See also
%    blk_diag_apply

function h_fb = radial_filter_to_fb_mat(h_fun, basis)
    if basis.type ~= fb_basis_type()
        error('Basis `basis` must be of Fourier-Bessel type.');
    end

    if basis.domain ~= 1
        error('Basis `basis` must be in the Fourier domain.');
    end

    if nargin(h_fun) ~= 1
        error('Function handle `h_fun` must take only one input.');
    end

    [k_vals, wts] = lgwt(basis.sz(1), 0, 0.5);

    h_vals = h_fun(k_vals);

    ind = 1;

    for ell = 0:basis.ell_max
        k_max = basis.k_max(ell+1);

        fb_vals = besselj(ell, 2*k_vals*basis.r0(1:k_max,ell+1)');
        fb_nrms = 1/sqrt(2)*abs(besselj(ell+1, basis.r0(1:k_max,ell+1)'))/2;
        fb_vals = bsxfun(@times, fb_vals, 1./fb_nrms);

        h_fb_vals = bsxfun(@times, h_vals, fb_vals);

        h_fb{ind} = fb_vals'*bsxfun(@times, k_vals.*wts, h_fb_vals);
        ind = ind+1;

        if ell > 0
            h_fb{ind} = h_fb{ind-1};
            ind = ind+1;
        end
    end
end
