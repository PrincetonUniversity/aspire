% FB_BASIS Construct a Fourier-Bessel basis object
%
% Usage
%    basis = fb_basis(sz, ell_max, domain);
%
% Input
%    sz: The size of the vectors for which to define the basis. Currently
%       only square images and cubic volumes are supported.
%    ell_max: The maximum order ell of the basis elements. If set to Inf,
%       the basis includes all ell such that the resulting basis vectors
%       are concentrated below the Nyquist frequency (default Inf).
%    domain: Specifies whether the decomposition should be in the spatial (0)
%       or frequency (1) domain (default 0).
%
% Output
%    basis: A Fourier-Bessel basis object corresponding to the parameters.
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
%       basis = fb_basis(sz, floor(sz(1)/2));
%       x = randn(sz);
%       v = basis.expand(x);
%
%    Likewise, we can map a random coefficient vector into the standard basis
%    using
%
%       v = randn(basis.count, 1);
%       x = basis.evaluate(v);
%
%    The adjoint of the evaluate function is obtained through basis.evaluate_t.
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

function basis = fb_basis(sz, ell_max, domain)
    if nargin < 2 || isempty(ell_max)
        ell_max = Inf;
    end

    if nargin < 3 || isempty(domain)
        domain = 0;
    end

    if numel(sz) ~= 3 && numel(sz) ~= 2
        error('Only two- or three-dimesional basis functions are supported.');
    end

    if ~all(sz == sz(1))
        error('Only cubic domains are supported.');
    end

    d = numel(sz);

    N = sz(1);

    % 2*sz(1) is an upper bound for ell_max, so use it here.
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

    basis.type = fb_basis_type();

    basis.sz = sz;

    basis.domain = domain;

    if d == 2
        basis.count = k_max(1) + sum(2*k_max(2:end));
    elseif d == 3
        basis.count = sum(k_max.*(2*[0:ell_max]'+1));
    end

    basis.ell_max = ell_max;
    basis.k_max = k_max;
    basis.r0 = r0;

    basis.indices = fb_indices(basis);

    basis.precomp = fb_precomp(basis);

    if isempty(basis.precomp)
        basis = rmfield(basis, 'precomp');
    end

    basis.evaluate = @(v)(fb_evaluate(v, basis));
    basis.evaluate_t = @(x)(fb_evaluate_t(x, basis));
    basis.expand = @(x)(fb_expand(x, basis));

    d = numel(basis.sz);

    basis.mat_evaluate = @(V)(mdim_mat_fun_conj(V, 1, d, basis.evaluate));
    basis.mat_expand = @(X)(mdim_mat_fun_conj(X, d, 1, basis.expand));
    basis.mat_evaluate_t = @(X)(mdim_mat_fun_conj(X, d, 1, basis.evaluate_t));
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

function [r_unique, ang_unique, r_idx, ang_idx, mask, stat_mask] = ...
    unique_coordinates_2d(N)

    mesh2d = mesh_2d(N);

    [mask, stat_mask] = positive_fourier_mask(N*ones(1, 2));
    mask = mdim_ifftshift(mask, 1:2);
    stat_mask = mdim_ifftshift(stat_mask, 1:2);

    mask = mask & mesh2d.r<=1;

    r = mesh2d.r(mask);
    phi = mesh2d.phi(mask);

    [r_unique, ~, r_idx] = unique(r);
    [ang_unique, ~, ang_idx] = unique(phi);
end

function [r_unique, ang_unique, r_idx, ang_idx, mask, stat_mask] = ...
    unique_coordinates_3d(N)

    mesh3d = mesh_3d(N);

    [mask, stat_mask] = positive_fourier_mask(N*ones(1, 3));
    mask = mdim_ifftshift(mask, 1:3);
    stat_mask = mdim_ifftshift(stat_mask, 1:3);

    mask = mask & mesh3d.r<=1;

    r = mesh3d.r(mask);
    theta = mesh3d.theta(mask);
    phi = mesh3d.phi(mask);

    [r_unique, ~, r_idx] = unique(r);
    [ang_unique, ~, ang_idx] = unique([theta phi], 'rows');
end

function nrm = basis_norm_2d(basis, ell, k)
    nrm = abs(besselj(ell+1, basis.r0(k,ell+1)))*sqrt(pi/2)* ...
        sqrt(prod(basis.sz/2));

    if ell == 0
        nrm = nrm*sqrt(2);
    end
end

function nrm = basis_norm_3d(basis, ell, k)
    nrm = abs(sph_bessel(ell+1, basis.r0(k,ell+1)))/sqrt(2)* ...
        sqrt(prod(basis.sz/2));
end

function precomp = fb_precomp(basis)
    d = numel(basis.sz);

    if d == 2
        [r_unique, ang_unique, r_idx, ang_idx, mask] = ...
            unique_coordinates_2d(basis.sz(1));

        ind_radial = 1;
        ind_ang = 1;

        radial = zeros(size(r_unique, 1), sum(basis.k_max));
        ang = zeros(size(ang_unique, 1), 1+2*(basis.ell_max-1));

        for ell = 0:basis.ell_max
            for k = 1:basis.k_max(ell+1)
                radial(:,ind_radial) = ...
                    besselj(ell, basis.r0(k,ell+1)*r_unique);

                ind_radial = ind_radial+1;
            end

            if ell == 0
                sgns = 1;
            else
                sgns = [1 -1];
            end

            for sgn = sgns
                if sgn == 1
                    ang(:,ind_ang) = cos(ell*ang_unique);
                else
                    ang(:,ind_ang) = sin(ell*ang_unique);
                end

                ind_ang = ind_ang+1;
            end
        end
    elseif d == 3
        [r_unique, ang_unique, r_idx, ang_idx, mask] = ...
            unique_coordinates_3d(basis.sz(1));

        ind_radial = 1;
        ind_ang = 1;

        radial = zeros(size(r_unique, 1), sum(basis.k_max));
        ang = zeros(size(ang_unique, 1), sum(2:[0:basis.ell_max]+1));

        for ell = 0:basis.ell_max
            for k = 1:basis.k_max(ell+1)
                radial(:,ind_radial) = ...
                    sph_bessel(ell, basis.r0(k,ell+1)*r_unique);

                ind_radial = ind_radial+1;
            end

            for m = -ell:ell
                ang(:,ind_ang) = real_sph_harmonic(ell, m, ...
                    ang_unique(:,1), ang_unique(:,2));

                ind_ang = ind_ang+1;
            end
        end
    end

    precomp = struct();
    precomp.radial = radial;
    precomp.ang = ang;
end

function x = fb_evaluate(v, basis)
    basis_check_evaluate(v, basis);

    d = numel(basis.sz);

    if d == 2
        x = fb_evaluate_2d(v, basis);
    elseif d == 3
        x = fb_evaluate_3d(v, basis);
    end
end

function x = fb_evaluate_2d(v, basis)
    [v, sz_roll] = unroll_dim(v, 2);

    [r_unique, ang_unique, r_idx, ang_idx, mask, stat_mask] = ...
        unique_coordinates_2d(basis.sz(1));

    is_precomp = isfield(basis, 'precomp');

    is_fourier = (isfield(basis, 'domain') && basis.domain == 1);

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    x_even = zeros([prod(basis.sz) size(v, 2)], class(v));
    x_odd = zeros([prod(basis.sz) size(v, 2)], class(v));
    for ell = 0:basis.ell_max
        k_max = basis.k_max(ell+1);

        idx_radial = ind_radial + [0:k_max-1];

        nrms = zeros(k_max, 1);
        for k = 1:k_max
            nrms(k) = basis_norm_2d(basis, ell, k);
        end

        if ~is_precomp
            radial = zeros(size(r_unique, 1), k_max);
            for k = 1:k_max
                radial(:,k) = besselj(ell, basis.r0(k,ell+1)*r_unique);
            end
        else
            radial = basis.precomp.radial(:,idx_radial);
        end

        radial = bsxfun(@times, radial, 1./nrms');

        if ell == 0
            sgns = 1;
        else
            sgns = [1 -1];
        end

        for sgn = sgns
            if ~is_precomp
                if sgn == 1
                    ang = cos(ell*ang_unique);
                else
                    ang = sin(ell*ang_unique);
                end
            else
                ang = basis.precomp.ang(:,ind_ang);
            end

            ang_radial = bsxfun(@times, ang(ang_idx), radial(r_idx,:));

            idx = ind + [0:k_max-1];

            if mod(ell, 2) == 0
                x_even(mask,:) = x_even(mask,:) + ang_radial*v(idx,:);
            else
                x_odd(mask,:) = x_odd(mask,:) + ang_radial*v(idx,:);
            end

            ind = ind + numel(idx);

            ind_ang = ind_ang + 1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    x_even(stat_mask,:) = 1/2*x_even(stat_mask,:);

    x_even = reshape(x_even, [basis.sz size(x_even, 2)]);
    x_odd = reshape(x_odd, [basis.sz size(x_odd, 2)]);

    x_even = x_even + fourier_flip(x_even, [1 2]);
    x_odd = x_odd - fourier_flip(x_odd, [1 2]);

    if is_fourier
        x_f = x_even + 1i*x_odd;
        x = sqrt(prod(basis.sz))*real(icfft2(x_f));
    else
        x = x_even + x_odd;
    end

    x = roll_dim(x, sz_roll);
end

function x = fb_evaluate_3d(v, basis)
    [v, sz_roll] = unroll_dim(v, 2);

    [r_unique, ang_unique, r_idx, ang_idx, mask, stat_mask] = ...
        unique_coordinates_3d(basis.sz(1));

    is_precomp = isfield(basis, 'precomp');

    is_fourier = (isfield(basis, 'domain') && basis.domain == 1);

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    x_even = zeros([prod(basis.sz) size(v, 2)], class(v));
    x_odd = zeros([prod(basis.sz) size(v, 2)], class(v));
    for ell = 0:basis.ell_max
        idx_radial = ind_radial + [0:basis.k_max(ell+1)-1];

        nrms = zeros(numel(idx_radial), 1);
        for k = 1:numel(idx_radial)
            nrms(k) = basis_norm_3d(basis, ell, k);
        end

        if ~is_precomp
            radial = zeros(size(r_unique, 1), numel(idx_radial));
            for k = 1:numel(idx_radial)
                radial(:,k) = sph_bessel(ell, basis.r0(k,ell+1)*r_unique);
            end
        else
            radial = basis.precomp.radial(:,idx_radial);
        end

        radial = bsxfun(@times, radial, 1./nrms');

        for m = -ell:ell
            if ~is_precomp
                ang = real_sph_harmonic(ell, m, ang_unique(:,1), ang_unique(:,2));
            else
                ang = basis.precomp.ang(:,ind_ang);
            end

            ang_radial = bsxfun(@times, ang(ang_idx), radial(r_idx,:));

            idx = ind + [0:numel(idx_radial)-1];

            if mod(ell, 2) == 0
                x_even(mask,:) = x_even(mask,:) + ang_radial*v(idx,:);
            else
                x_odd(mask,:) = x_odd(mask,:) + ang_radial*v(idx,:);
            end

            ind = ind + numel(idx);
            ind_ang = ind_ang+1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    x_even(stat_mask,:) = 1/2*x_even(stat_mask,:);

    x_even = reshape(x_even, [basis.sz size(x_even, 2)]);
    x_odd = reshape(x_odd, [basis.sz size(x_odd, 2)]);

    x_even = x_even + fourier_flip(x_even, 1:3);
    x_odd = x_odd - fourier_flip(x_odd, 1:3);

    if is_fourier
        x_f = x_even + 1i*x_odd;
        x = sqrt(prod(basis.sz))*real(icfft3(x_f));
    else
        x = x_even + x_odd;
    end

    x = roll_dim(x, sz_roll);
end

function v = fb_evaluate_t(x, basis)
    basis_check_expand(x, basis);

    d = numel(basis.sz);

    if d == 2
        v = fb_evaluate_t_2d(x, basis);
    elseif d == 3
        v = fb_evaluate_t_3d(x, basis);
    end
end

function v = fb_evaluate_t_2d(x, basis)
    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    [r_unique, ang_unique, r_idx, ang_idx, mask, stat_mask] = ...
        unique_coordinates_2d(basis.sz(1));

    is_fourier = (isfield(basis, 'domain') && basis.domain == 1);

    if is_fourier
        x_f = 1/sqrt(prod(basis.sz))*cfft2(x);

        x_even = real(x_f)*2;
        x_odd = imag(x_f)*2;
    else
        x_even = (x + fourier_flip(x, [1 2]));
        x_odd = (x - fourier_flip(x, [1 2]));
    end

    x_even = reshape(x_even, [prod(basis.sz) size(x_even, numel(basis.sz)+1)]);
    x_odd = reshape(x_odd, [prod(basis.sz) size(x_odd, numel(basis.sz)+1)]);

    x_even(stat_mask,:) = 1/2*x_even(stat_mask,:);

    is_precomp = isfield(basis, 'precomp');

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    v = zeros([basis.count size(x_even, 2)], class(x));

    for ell = 0:basis.ell_max
        idx_radial = ind_radial + [0:basis.k_max(ell+1)-1];

        nrms = zeros(numel(idx_radial), 1);
        for k = 1:numel(idx_radial)
            nrms(k) = basis_norm_2d(basis, ell, k);
        end

        if ~is_precomp
            radial = zeros(size(r_unique, 1), numel(idx_radial));
            for k = 1:numel(idx_radial)
                radial(:,k) = besselj(ell, basis.r0(k,ell+1)*r_unique);
            end
        else
            radial = basis.precomp.radial(:,idx_radial);
        end

        radial = bsxfun(@times, radial, 1./nrms');

        if ell == 0
            sgns = 1;
        else
            sgns = [1 -1];
        end

        for sgn = sgns
            if ~is_precomp
                if sgn == 1
                    ang = cos(ell*ang_unique);
                else
                    ang = sin(ell*ang_unique);
                end
            else
                ang = basis.precomp.ang(:,ind_ang);
            end

            ang_radial = bsxfun(@times, ang(ang_idx), radial(r_idx,:));

            idx = ind + [0:numel(idx_radial)-1];

            if mod(ell, 2) == 0
                v(idx,:) = ang_radial'*x_even(mask,:);
            else
                v(idx,:) = ang_radial'*x_odd(mask,:);
            end

            ind = ind + numel(idx);

            ind_ang = ind_ang+1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    v = roll_dim(v, sz_roll);
end

function v = fb_evaluate_t_3d(x, basis)
    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    [r_unique, ang_unique, r_idx, ang_idx, mask, stat_mask] = ...
        unique_coordinates_3d(basis.sz(1));

    is_fourier = (isfield(basis, 'domain') && basis.domain == 1);

    if is_fourier
        x_f = 1/sqrt(prod(basis.sz))*cfft3(x);

        x_even = 2*real(x_f);
        x_odd = 2*imag(x_f);
    else
        x_even = x + fourier_flip(x, 1:3);
        x_odd = x - fourier_flip(x, 1:3);
    end

    x_even = reshape(x_even, [prod(basis.sz) size(x_even, numel(basis.sz)+1)]);
    x_odd = reshape(x_odd, [prod(basis.sz) size(x_odd, numel(basis.sz)+1)]);

    x_even(stat_mask,:) = 1/2*x_even(stat_mask,:);

    is_precomp = isfield(basis, 'precomp');

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    v = zeros([basis.count size(x_even, 2)], class(x));
    for ell = 0:basis.ell_max
        idx_radial = ind_radial + [0:basis.k_max(ell+1)-1];

        nrms = zeros(numel(idx_radial), 1);
        for k = 1:numel(idx_radial)
            nrms(k) = basis_norm_3d(basis, ell, k);
        end

        if ~is_precomp
            radial = zeros(size(r_unique, 1), numel(idx_radial));
            for k = 1:numel(idx_radial)
                radial(:,k) = sph_bessel(ell, basis.r0(k,ell+1)*r_unique);
            end
        else
            radial = basis.precomp.radial(:,idx_radial);
        end

        radial = bsxfun(@times, radial, 1./nrms');

        for m = -ell:ell
            if ~is_precomp
                ang = real_sph_harmonic(ell, m, ang_unique(:,1), ang_unique(:,2));
            else
                ang = basis.precomp.ang(:,ind_ang);
            end

            ang_radial = bsxfun(@times, ang(ang_idx), radial(r_idx,:));

            idx = ind + [0:numel(idx_radial)-1];

            if mod(ell, 2) == 0
                v(idx,:) = ang_radial'*x_even(mask,:);
            else
                v(idx,:) = ang_radial'*x_odd(mask,:);
            end

            ind = ind + numel(idx);
            ind_ang = ind_ang+1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    v = roll_dim(v, sz_roll);
end

function v = fb_expand(x, basis)
    basis_check_expand(x, basis);

    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    b = basis.evaluate_t(x);

    A = @(v)(basis.evaluate_t(basis.evaluate(v)));

    % TODO: Check that this tolerance make sense for multiple columns in x

    cg_opt.max_iter = Inf;
    cg_opt.rel_tolerance = 10*eps(class(x));
    cg_opt.verbose = false;

    v = conj_grad(A, b, cg_opt);

    v = roll_dim(v, sz_roll);
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
