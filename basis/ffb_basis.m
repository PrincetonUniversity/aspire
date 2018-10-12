% FFB_BASIS Construct a fast Fourier-Bessel basis object
%
% Usage
%    basis = ffb_basis(sz, ell_max, domain);
%
% Input
%    sz: The size of the vectors for which to define the basis. Currently
%       only square images are supported.
%    ell_max: The maximum order ell of the basis elements. If set to Inf,
%       the basis includes all ell such that the resulting basis vectors
%       are concentrated below the Nyquist frequency (default Inf).
%    domain: Specifies whether the decomposition should be in the spatial (0)
%       or frequency (1) domain. Currently, only frequency-domain basis
%       vectors are supported (default 1).
%
% Output
%    basis: A fast Fourier-Bessel basis object corresponding to the parameters.
%
% Description
%    The returned basis object can be used to project a vector x in the
%    standard basis onto the basis or to map a coefficient vector v into
%    the standard basis. This is achieved using the basis.expand and
%    basis.evaluate functions, respectively. The fields basis.sz denotes
%    the size of the vectors in the standard basis while basis.count denotes
%    the number of vectors in the basis.
%
%    For example, we can generate a random vector in the standard basis and
%    project it onto the Fourier-Bessel basis through
%
%       sz = 16*ones(1, 2);
%       basis = ffb_basis(sz, floor(sz(1)/2));
%       x = randn(sz);
%       v = basis.expand(x);
%
%    Likewise, we can map a random coefficient vector into the standard basis
%    using
%
%       v = randn(basis.count, 1);
%       x = basis.evaluate(v);
%
%    The adjoint of the evaluate function is obtained through
%    basis.evaluate_t.
%
%    This particular basis is a Fourier-Bessel basis in two or three
%    dimensions. It is a separable basis consisting of a radial part multiplied
%    by an angular part.
%
%    The radial part is a Bessel function in 2D or a spherical Bessel function
%    in 3D. These are dilated to obtain a given number of zeros on the
%    interval [0, 1]. Specifically, they are of the form
%
%       f_{ell,k}(r) = Z_ell(R_{ell,k} r),
%
%    where Z_ell is the ellth-order (spherical) Bessel function and R_{ell,k}
%    is its kth zero.
%
%    The angular part is given by a sinusoid in 2D or a real spherical harmonic
%    function in 3D. Specifically, for ell = 0, the angular part is constant,
%    while for ell > 0, it is either a sine or a cosine.
%
%    All of the basis vectors are normalized so that their 2-norm, as continuous
%    functions, is equal to 1. Additionally, the basis functions are orthogonal,
%    again in the continuous case. Due to the discrete sampling, these
%    properties only hold asymptotically as max(sz) goes to infinity.
%
%    The parameter ell_max determines the highest frequency of the angular
%    sinusoids in 2D and the highest degree of the real spherical harmonics in
%    3D. For each ell, the maximum number of zeros in the radial part f_{ell,k}
%    is determined by ensuring that the Fourier transform of the entire basis
%    function is mostly concentrated within the Nyquist disk or ball. This
%    maximum number is denoted by k_max(ell) and is stored in the basis object.
%
%    The coefficients of the basis are ordered by increasing ell, from 0 to
%    ell_max. For each ell, the basis functions are ordered by the radial index
%    k from 1 to k_max(ell). In 2D, for each radial index, we then have the
%    cosine part followed by the sine part (except for ell = 0, in which case we
%    only have the cosine part). In 3D, we instead have the coefficients
%    arranged by order m, ranging from -ell to +ell.
%
% Note
%    The underlying continuous basis functions for ffb_basis is the same as
%    for fb_basis, but they differ in their discretizations. In particular,
%    the discretization of ffb_basis allows for the use of NUFFTs, which
%    significantly speeds up the evaluation and adjoint mappings.

function basis = ffb_basis(sz, ell_max, domain)
    if nargin < 2 || isempty(ell_max)
        ell_max = Inf;
    end

    if nargin < 3 || isempty(domain)
        domain = 1;
    end

    if ~all(sz == sz(1))
        error('Only square/cubic domains are supported.');
    end

    if domain == 0
        error('Only frequency domain is supported.');
    end

    d = numel(sz);

    N = sz(1);

    k_max = zeros(min(ell_max+1, 2*sz(1)+1), 1);

    r0 = cell(1, numel(k_max));

    for ell = 0:numel(k_max)-1
        [k_max(ell+1), r0{ell+1}] = num_besselj_zeros(ell+(d-2)/2, N*pi/2);

        if k_max(ell+1) == 0
            ell_max = ell-1;

            k_max = k_max(1:ell_max+1);
            r0 = r0(1:ell_max+1);
            break;
        end
    end

    for ell = 0:ell_max
        r0{ell+1} = cat(1, r0{ell+1}', NaN(max(k_max)-k_max(ell+1), 1));
    end

    r0 = cell2mat(r0);

    basis = struct();

    basis.type = ffb_basis_type();

    basis.sz = sz;

    basis.domain = domain;

    basis.R = N/2;
    basis.c = 0.5;

    if d == 2
        basis.count = k_max(1) + sum(2*k_max(2:end));
    elseif d == 3
        basis.count = sum(k_max.*(2*[0:ell_max]'+1));
    end

    basis.ell_max = ell_max;
    basis.k_max = k_max;
    basis.r0 = r0;

    basis.indices = fb_indices(basis);

    basis.precomp = ffb_precomp(basis);

    basis.evaluate = @(v)(ffb_evaluate(v, basis));
    basis.evaluate_t = @(x)(ffb_evaluate_t(x, basis));

    basis.expand = @(x)(ffb_expand(x, basis));

    d = numel(basis.sz);

    basis.mat_evaluate = @(V)(mdim_mat_fun_conj(V, 1, d, basis.evaluate));
    basis.mat_evaluate_t = @(X)(mdim_mat_fun_conj(X, d, 1, basis.evaluate_t));

    basis.mat_expand = @(X)(mdim_mat_fun_conj(X, d, 1, basis.expand));
end

function [k, r0] = num_besselj_zeros(ell, r)
    k = 4;

    r0 = besselj_zeros(ell, k);
    while all(r0 < r)
        k = 2*k;

        r0 = besselj_zeros(ell, k);
    end

    r0 = r0(r0 < r);

    k = numel(r0);
end

function nrm = basis_norm_2d(basis, ell, k)
    nrm = abs(besselj(ell+1, basis.r0(k,ell+1)))*sqrt(pi/2)* ...
        sqrt(prod(basis.sz/2));

    if ell == 0
        nrm = nrm*sqrt(2);
    end
end

function indices = fb_indices(basis)
    indices = struct();

    d = numel(basis.sz);

    if d == 2
        indices.ells = zeros(basis.count, 1);
        indices.ks = zeros(basis.count, 1);
        indices.sgns = zeros(basis.count, 1);

        ind = 1;

        for ell = 0:basis.ell_max
            if ell == 0
                sgns = 1;
            else
                sgns = [1 -1];
            end

            ks = 1:basis.k_max(ell+1);

            for sgn = sgns
                rng = ind:ind+numel(ks)-1;

                indices.ells(rng) = ell;
                indices.ks(rng) = ks;
                indices.sgns(rng) = sgn;

                ind = ind + numel(rng);
            end
        end
    elseif d == 3
        indices.ells = zeros(basis.count, 1);
        indices.ms = zeros(basis.count, 1);
        indices.ks = zeros(basis.count, 1);

        ind = 1;

        for ell = 0:basis.ell_max
            ks = 1:basis.k_max(ell+1);
            ms = -ell:ell;

            for m = -ell:ell
                rng = ind:ind+numel(ks)-1;

                indices.ells(rng) = ell;
                indices.ms(rng) = m;
                indices.ks(rng) = ks;

                ind = ind+numel(ks);
            end
        end
    end
end

function precomp = ffb_precomp(basis)
    d = numel(basis.sz);

    if d == 2
        precomp = ffb_precomp_2d(basis);
    elseif d == 3
        precomp = ffb_precomp_3d(basis);
    end
end

function precomp = ffb_precomp_2d(basis)
    n_r = ceil(4*basis.c*basis.R);

    [precomp.r, precomp.w] = lgwt(n_r, 0, basis.c);

    % Second dimension below is not basis.count because we're not counting
    % signs, since radial part is the same for both cos and sin.
    precomp.radial = zeros(n_r, sum(basis.k_max));

    ind = 1;

    for ell = 0:basis.ell_max
        idx = ind + [0:basis.k_max(ell+1)-1];

        besselj_zeros = basis.r0(1:basis.k_max(ell+1), ell+1)';
        radial_ell = besselj(ell, 2*precomp.r*besselj_zeros);

        % NOTE: We need to remove the factor due to the discretization here
        % since it is already included in our quadrature weights.
        nrms = 1/(sqrt(prod(basis.sz)))*basis_norm_2d(basis, ell, [1:basis.k_max(ell+1)])';

        radial_ell = bsxfun(@times, radial_ell, 1./nrms);

        precomp.radial(:,idx) = radial_ell;

        ind = ind + numel(idx);
    end

    n_theta = ceil(16*basis.c*basis.R);
    n_theta = (n_theta + mod(n_theta, 2))/2;

    % Only calculate "positive" frequencies in one half-plane.
    freqs_x = precomp.r*cos([0:n_theta-1]*2*pi/(2*n_theta));
    freqs_y = precomp.r*sin([0:n_theta-1]*2*pi/(2*n_theta));

    precomp.freqs = cat(1, permute(freqs_x, [3 1 2]), permute(freqs_y, [3 1 2]));
end

function precomp = ffb_precomp_3d(basis)
    precomp = struct();

    n_r = basis.sz(1);
    n_theta = 2*basis.sz(1);
    n_phi = basis.sz(1);

    [r, wt_r] = lgwt(n_r, 0, 1/2);
    [z, wt_z] = lgwt(n_phi, -1, 1);
    phi = acos(z);
    wt_phi = wt_z;
    theta = 2*pi*[0:n_theta-1]'/(2*n_theta);

    precomp.radial_wtd = zeros([n_r max(basis.k_max) basis.ell_max+1]);
    for ell = 0:basis.ell_max
        k_max_ell = basis.k_max(ell+1);

        radial_ell = sph_bessel(ell, 2*r*basis.r0(1:k_max_ell,ell+1)');
        nrm = abs(sph_bessel(ell+1, basis.r0(1:k_max_ell,ell+1)')/4);
        radial_ell = bsxfun(@times, radial_ell, 1./nrm);

        radial_ell_wtd = bsxfun(@times, r.^2.*wt_r, radial_ell);

        precomp.radial_wtd(:,1:k_max_ell,ell+1) = radial_ell_wtd;
    end

    precomp.angular_phi_wtd_even = {};
    precomp.angular_phi_wtd_odd = {};

    for m = 0:basis.ell_max
        n_even_ell = floor((basis.ell_max-m)/2)+1 ...
            -mod(basis.ell_max, 2)*mod(m, 2);
        n_odd_ell = basis.ell_max-m+1-n_even_ell;

        angular_phi_wtd_m_even = zeros(n_phi, n_even_ell);
        angular_phi_wtd_m_odd = zeros(n_phi, n_odd_ell);

        ind_even = 1;
        ind_odd = 1;

        for ell = m:basis.ell_max
            angular_phi_m_ell = norm_assoc_legendre(ell, m, cos(phi));

            nrm_inv = sqrt(1/(2*pi));
            angular_phi_m_ell = angular_phi_m_ell*nrm_inv;

            angular_phi_wtd_m_ell = wt_phi.*angular_phi_m_ell;

            if mod(ell, 2) == 0
                angular_phi_wtd_m_even(:,ind_even) = ...
                    angular_phi_wtd_m_ell;
                ind_even = ind_even+1;
            else
                angular_phi_wtd_m_odd(:,ind_odd) = ...
                    angular_phi_wtd_m_ell;
                ind_odd = ind_odd+1;
            end
        end

        precomp.angular_phi_wtd_even{m+1} = angular_phi_wtd_m_even;
        precomp.angular_phi_wtd_odd{m+1} = angular_phi_wtd_m_odd;
    end

    angular_theta = zeros(n_theta, 2*basis.ell_max+1);

    angular_theta(:,1:basis.ell_max) = ...
        sqrt(2)*sin(theta*[basis.ell_max:-1:1]);
    angular_theta(:,basis.ell_max+1) = ones(n_theta, 1);
    angular_theta(:,basis.ell_max+2:2*basis.ell_max+1) = ...
        sqrt(2)*cos(theta*[1:basis.ell_max]);

    angular_theta_wtd = 2*pi/n_theta*angular_theta;

    precomp.angular_theta_wtd = angular_theta_wtd;

    [theta_grid, phi_grid, r_grid] = ndgrid(theta, phi, r);

    fourier_x = r_grid .* cos(theta_grid) .* sin(phi_grid);
    fourier_y = r_grid .* sin(theta_grid) .* sin(phi_grid);
    fourier_z = r_grid .* cos(phi_grid);

    precomp.fourier_pts = 2*pi*[fourier_x(:) fourier_y(:) fourier_z(:)]';
end

function x = ffb_evaluate(v, basis)
    basis_check_evaluate(v, basis);

    d = numel(basis.sz);

    if d == 2
        x = ffb_evaluate_2d(v, basis);
    elseif d == 3
        x = ffb_evaluate_3d(v, basis);
    end
end

function x = ffb_evaluate_2d(v, basis)
    n_theta = size(basis.precomp.freqs, 3);
    n_r = size(basis.precomp.freqs, 2);
    n_data = size(v, 2);

    % TODO: Rename. This is not actually the polar FT.
    pf = zeros([n_r 2*n_theta n_data], class(v));

    mask = (basis.indices.ells == 0);

    ind = 1;

    idx = ind + [0:basis.k_max(1)-1];

    pf(:,1,:) = basis.precomp.radial(:,idx)*v(mask,:);

    ind = ind + numel(idx);

    ind_pos = ind;

    for ell = 1:basis.ell_max
        idx = ind + [0:basis.k_max(ell+1)-1];

        idx_pos = ind_pos + [0:basis.k_max(ell+1)-1];
        idx_neg = idx_pos + basis.k_max(ell+1);

        v_ell = (v(idx_pos,:) - 1i*v(idx_neg,:))/2;

        if mod(ell, 2) == 1
            v_ell = 1i*v_ell;
        end

        pf_ell = basis.precomp.radial(:,idx)*v_ell;
        pf(:,ell+1,:) = pf_ell;

        if mod(ell, 2) == 0
            pf(:,end-ell+1,:) = conj(pf_ell);
        else
            pf(:,end-ell+1,:) = -conj(pf_ell);
        end

        ind = ind + numel(idx);

        ind_pos = ind_pos + 2*basis.k_max(ell+1);
    end

    pf = 2*pi*ifft(pf, [], 2);

    % Only need "positive" frequencies.
    pf = pf(:,1:end/2,:);

    pf = bsxfun(@times, pf, basis.precomp.w.*basis.precomp.r);

    pf = reshape(pf, [n_r*n_theta n_data]);
    freqs = reshape(basis.precomp.freqs, [2 n_r*n_theta]);

    x = 2*real(anufft2(pf, 2*pi*freqs, basis.sz));
end

function x = ffb_evaluate_3d(v, basis)
    n_data = size(v, 2);

    n_r = size(basis.precomp.radial_wtd, 1);
    n_phi = size(basis.precomp.angular_phi_wtd_even{1}, 1);
    n_theta = size(basis.precomp.angular_theta_wtd, 1);

    u_even = zeros( ...
        [n_r 2*basis.ell_max+1 n_data floor(basis.ell_max/2)+1], class(v));
    u_odd = zeros( ...
        [n_r 2*basis.ell_max+1 n_data ceil(basis.ell_max/2)], class(v));

    for ell = 0:basis.ell_max
        k_max_ell = basis.k_max(ell+1);

        radial_wtd = basis.precomp.radial_wtd(:,1:k_max_ell,ell+1);

        % TODO: Fix this to avoid lookup each time.
        ind = (basis.indices.ells == ell);

        v_ell = reshape(v(ind,:), [k_max_ell (2*ell+1)*n_data]);

        v_ell = radial_wtd*v_ell;

        v_ell = reshape(v_ell, [n_r 2*ell+1 n_data]);

        if mod(ell, 2) == 0
            u_even(:,[-ell:ell]+basis.ell_max+1,:,ell/2+1) = v_ell;
        else
            u_odd(:,[-ell:ell]+basis.ell_max+1,:,(ell-1)/2+1) = v_ell;
        end
    end

    u_even = permute(u_even, [4 1 2 3]);
    u_odd = permute(u_odd, [4 1 2 3]);

    w_even = zeros([n_phi n_r n_data 2*basis.ell_max+1], class(v));
    w_odd = zeros([n_phi n_r n_data 2*basis.ell_max+1], class(v));

    for m = 0:basis.ell_max
        angular_phi_wtd_m_even = basis.precomp.angular_phi_wtd_even{m+1};
        angular_phi_wtd_m_odd = basis.precomp.angular_phi_wtd_odd{m+1};

        n_even_ell = size(angular_phi_wtd_m_even, 2);
        n_odd_ell = size(angular_phi_wtd_m_odd, 2);

        if m == 0
            sgns = 1;
        else
            sgns = [1 -1];
        end

        for sgn = sgns
            u_m_even = u_even(end-n_even_ell+1:end,:,basis.ell_max+1+sgn*m,:);
            u_m_odd = u_odd(end-n_odd_ell+1:end,:,basis.ell_max+1+sgn*m,:);

            u_m_even = reshape(u_m_even, [n_even_ell n_r*n_data]);
            u_m_odd = reshape(u_m_odd, [n_odd_ell n_r*n_data]);

            w_m_even = angular_phi_wtd_m_even*u_m_even;
            w_m_odd = angular_phi_wtd_m_odd*u_m_odd;

            w_m_even = reshape(w_m_even, [n_phi n_r n_data]);
            w_m_odd = reshape(w_m_odd, [n_phi n_r n_data]);

            w_even(:,:,:,basis.ell_max+1+sgn*m) = w_m_even;
            w_odd(:,:,:,basis.ell_max+1+sgn*m) = w_m_odd;
        end
    end

    w_even = permute(w_even, [4 1 2 3]);
    w_odd = permute(w_odd, [4 1 2 3]);

    u_even = w_even;
    u_odd = w_odd;

    u_even = reshape(u_even, [2*basis.ell_max+1 n_phi*n_r*n_data]);
    u_odd = reshape(u_odd, [2*basis.ell_max+1 n_phi*n_r*n_data]);

    w_even = basis.precomp.angular_theta_wtd*u_even;
    w_odd = basis.precomp.angular_theta_wtd*u_odd;

    pf = w_even + 1i*w_odd;
    pf = reshape(pf, [n_theta*n_phi*n_r n_data]);

    x = real(anufft3(pf, basis.precomp.fourier_pts, basis.sz));

    x = cast(x, class(v));
end

function v = ffb_evaluate_t(x, basis)
    basis_check_expand(x, basis);

    d = numel(basis.sz);

    if d == 2
        v = ffb_evaluate_t_2d(x, basis);
    elseif d == 3
        v = ffb_evaluate_t_3d(x, basis);
    end
end

function v = ffb_evaluate_t_2d(x, basis)
    n_theta = size(basis.precomp.freqs, 3);
    n_r = size(basis.precomp.freqs, 2);
    n_data = size(x, 3);

    freqs = reshape(basis.precomp.freqs, [2 n_r*n_theta]);
    pf = nufft2(x, 2*pi*freqs);
    pf = reshape(pf, [n_r n_theta n_data]);

    % Recover "negative" frequencies from "positive" half plane.
    pf = cat(2, pf, conj(pf));

    pf = bsxfun(@times, pf, basis.precomp.w.*basis.precomp.r);

    % TODO: Rename. This isn't actually the polar FT.
    pf = 2*pi/(2*n_theta)*fft(pf, [], 2);

    % This only makes it easier to slice the array later.
    pf = permute(pf, [1 3 2]);

    v = zeros([basis.count n_data], class(x));

    ind = 1;

    idx = ind + [0:basis.k_max(1)-1];

    mask = (basis.indices.ells == 0);

    v(mask,:) = basis.precomp.radial(:,idx)'*real(pf(:,:,1));

    ind = ind + numel(idx);

    ind_pos = ind;

    for ell = 1:basis.ell_max
        idx = ind + [0:basis.k_max(ell+1)-1];

        idx_pos = ind_pos + [0:basis.k_max(ell+1)-1];
        idx_neg = idx_pos + basis.k_max(ell+1);

        v_ell = basis.precomp.radial(:,idx)'*pf(:,:,ell+1);

        if mod(ell, 2) == 0
            v_pos = real(v_ell);
            v_neg = -imag(v_ell);
        else
            v_pos = imag(v_ell);
            v_neg = real(v_ell);
        end

        v(idx_pos,:) = v_pos;
        v(idx_neg,:) = v_neg;

        ind = ind + numel(idx);

        ind_pos = ind_pos + 2*basis.k_max(ell+1);
    end
end

function v = ffb_evaluate_t_3d(x, basis)
    n_data = size(x, 4);

    n_r = size(basis.precomp.radial_wtd, 1);
    n_phi = size(basis.precomp.angular_phi_wtd_even{1}, 1);
    n_theta = size(basis.precomp.angular_theta_wtd, 1);

    pf = nufft3(x, basis.precomp.fourier_pts);

    pf = cast(pf, class(x));

    pf = reshape(pf, [n_theta n_phi*n_r*n_data]);

    u_even = basis.precomp.angular_theta_wtd'*real(pf);
    u_odd = basis.precomp.angular_theta_wtd'*imag(pf);

    u_even = reshape(u_even, [2*basis.ell_max+1 n_phi n_r n_data]);
    u_odd = reshape(u_odd, [2*basis.ell_max+1 n_phi n_r n_data]);

    u_even = permute(u_even, [2 3 4 1]);
    u_odd = permute(u_odd, [2 3 4 1]);

    w_even = zeros( ...
        [floor(basis.ell_max/2)+1 n_r 2*basis.ell_max+1 n_data], class(x));
    w_odd = zeros( ...
        [ceil(basis.ell_max/2) n_r 2*basis.ell_max+1 n_data], class(x));

    for m = 0:basis.ell_max
        angular_phi_wtd_m_even = basis.precomp.angular_phi_wtd_even{m+1};
        angular_phi_wtd_m_odd = basis.precomp.angular_phi_wtd_odd{m+1};

        n_even_ell = size(angular_phi_wtd_m_even, 2);
        n_odd_ell = size(angular_phi_wtd_m_odd, 2);

        if m == 0
            sgns = 1;
        else
            sgns = [1 -1];
        end

        for sgn = sgns
            u_m_even = u_even(:,:,:,basis.ell_max+1+sgn*m);
            u_m_odd = u_odd(:,:,:,basis.ell_max+1+sgn*m);

            u_m_even = reshape(u_m_even, [n_phi n_r*n_data]);
            u_m_odd = reshape(u_m_odd, [n_phi n_r*n_data]);

            w_m_even = angular_phi_wtd_m_even'*u_m_even;
            w_m_odd = angular_phi_wtd_m_odd'*u_m_odd;

            w_m_even = reshape(w_m_even, [n_even_ell n_r n_data]);
            w_m_odd = reshape(w_m_odd, [n_odd_ell n_r n_data]);

            w_even(end-n_even_ell+1:end,:,basis.ell_max+1+sgn*m,:) = w_m_even;
            w_odd(end-n_odd_ell+1:end,:,basis.ell_max+1+sgn*m,:) = w_m_odd;
        end
    end

    w_even = permute(w_even, [2 3 4 1]);
    w_odd = permute(w_odd, [2 3 4 1]);

    v = zeros([basis.count n_data], class(x));

    for ell = 0:basis.ell_max
        k_max_ell = basis.k_max(ell+1);

        radial_wtd = basis.precomp.radial_wtd(:,1:k_max_ell,ell+1);

        if mod(ell, 2) == 0
            v_ell = w_even(:,[-ell:ell]+basis.ell_max+1,:,ell/2+1);
        else
            v_ell = w_odd(:,[-ell:ell]+basis.ell_max+1,:,(ell-1)/2+1);
        end

        v_ell = reshape(v_ell, [n_r (2*ell+1)*n_data]);

        v_ell = radial_wtd'*v_ell;

        v_ell = reshape(v_ell, [k_max_ell*(2*ell+1) n_data]);

        % TODO: Fix this to avoid lookup each time.
        ind = (basis.indices.ells == ell);

        v(ind,:) = v_ell;
    end
end

function v = ffb_expand(x, basis)
    basis_check_expand(x, basis);

    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    b = basis.evaluate_t(x);

    A = @(v)(basis.evaluate_t(basis.evaluate(v)));

    % TODO: Check that this tolerance make sense for multiple columns in x

    cg_opt.max_iter = Inf;
    cg_opt.rel_tolerance = 10*eps(class(x));
    cg_opt.verbose = 0;

    v = conj_grad(A, b, cg_opt);

    v = roll_dim(v, sz_roll);
end
