% FOURIER_BESSEL_BASIS Construct a Fourier-Bessel basis object
%
% Usage
%    basis = fourier_bessel_basis(sz, ell_max);
%
% Input
%    sz: The size of the vectors for which to define the basis. Currently
%       only cubic volumes are supported.
%    ell_max: The maximum order ell of the basis elements.
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
%       basis = fourier_bessel_basis(sz, mask);
%       x = randn([sz 1]);
%       v = basis.expand(x);
%
%    Likewise, we can map a random coefficient vector into the standard basis
%    using
%
%       v = randn(basis.count, 1);
%       x = basis.evaluate(v);
%
%    The adjoints of the expand and evaluate functions are obtained through
%    basis.expand_t and basis.evaluate_t functions, respectively.

function basis = fourier_bessel_basis(sz, ell_max)
    if numel(sz) ~= 3
        error('Only three-dimesional basis functions are supported.');
    end

    if ~all(sz == sz(1))
        error('Only cubic domains are supported.');
    end

    N = sz(1);

    k_max = zeros(ell_max+1, 1);

    r0 = cell(1, ell_max+1);

    for ell = 0:ell_max
        [k_max(ell+1), r0{ell+1}] = num_besselj_zeros(ell+1/2, N*pi/2);
    end

    for ell = 0:ell_max
        r0{ell+1} = cat(1, r0{ell+1}', NaN(max(k_max)-k_max(ell+1), 1));
    end

    r0 = cell2mat(r0);

    basis = struct();

    basis.sz = sz;
    basis.count = sum(k_max.*(2*[0:ell_max]'+1));

    basis.ell_max = ell_max;
    basis.k_max = k_max;
    basis.r0 = r0;

    basis.precomp = fourier_bessel_precomp(basis);

    basis.evaluate = @(v)(fourier_bessel_evaluate(v, basis));
    basis.evaluate_t = @(x)(fourier_bessel_evaluate_t(x, basis));
    basis.expand = @(x)(fourier_bessel_expand(x, basis));
    basis.expand_t = @(v)(fourier_bessel_expand_t(v, basis));

    d = numel(basis.sz);

    basis.mat_evaluate = @(V)(mdim_mat_fun_conj(V, 1, d, basis.evaluate));
    basis.mat_expand = @(X)(mdim_mat_fun_conj(X, d, 1, basis.expand));
    basis.mat_evaluate_t = @(X)(mdim_mat_fun_conj(X, d, 1, basis.evaluate_t));
    basis.mat_expand_t = @(V)(mdim_mat_fun_conj(V, 1, d, basis.expand_t));
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

function [r_unique, ang_unique, r_idx, ang_idx, mask] = unique_coordinates(N)
    mesh3d = mesh_3d(N);

    mask = mesh3d.r<=1;

    r = mesh3d.r(mask);
    theta = mesh3d.theta(mask);
    phi = mesh3d.phi(mask);

    [r_unique, ~, r_idx] = unique(r);
    [ang_unique, ~, ang_idx] = unique([theta phi], 'rows');
end

function nrm = basis_norm(basis, ell, k)
    nrm = abs(sph_bessel(ell+1, basis.r0(k,ell+1)))/sqrt(2)* ...
        sqrt(prod(basis.sz/2));
end

function precomp = fourier_bessel_precomp(basis)
    [r_unique, ang_unique, r_idx, ang_idx, mask] = ...
        unique_coordinates(basis.sz(1));

    ind_radial = 1;
    ind_ang = 1;

    radial = zeros(size(r_unique, 1), sum(basis.k_max));
    ang = zeros(size(ang_unique, 1), sum(2:[0:basis.ell_max]+1));

    for ell = 0:basis.ell_max
        for k = 1:basis.k_max(ell+1)
            radial(:,ind_radial) = sph_bessel(ell, basis.r0(k,ell+1)*r_unique);

            ind_radial = ind_radial+1;
        end

        for m = -ell:ell
            ang(:,ind_ang) = real_sph_harmonic(ell, m, ...
                ang_unique(:,1), ang_unique(:,2));

            ind_ang = ind_ang+1;
        end
    end

    precomp = struct();
    precomp.radial = radial;
    precomp.ang = ang;
end

function x = fourier_bessel_evaluate(v, basis)
    basis_check_evaluate(v, basis);

    [v, sz_roll] = unroll_dim(v, 2);

    [r_unique, ang_unique, r_idx, ang_idx, mask] = ...
        unique_coordinates(basis.sz(1));

    is_precomp = isfield(basis, 'precomp');

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    x = zeros([prod(basis.sz) size(v, 2)], class(v));
    for ell = 0:basis.ell_max
        idx_radial = ind_radial + [0:basis.k_max(ell+1)-1];

        nrms = zeros(numel(idx_radial), 1);
        for k = 1:numel(idx_radial)
            nrms(k) = basis_norm(basis, ell, k);
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

            x(mask,:) = x(mask,:) + ang_radial*v(idx,:);

            ind = ind + numel(idx);
            ind_ang = ind_ang+1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    x = reshape(x, [basis.sz size(x, 2)]);

    x = roll_dim(x, sz_roll);
end

function v = fourier_bessel_evaluate_t(x, basis)
    basis_check_expand(x, basis);

    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    x = reshape(x, [prod(basis.sz) size(x, numel(basis.sz)+1)]);

    [r_unique, ang_unique, r_idx, ang_idx, mask] = ...
        unique_coordinates(basis.sz(1));

    is_precomp = isfield(basis, 'precomp');

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    v = zeros([basis.count size(x, 2)], class(x));
    for ell = 0:basis.ell_max
        idx_radial = ind_radial + [0:basis.k_max(ell+1)-1];

        nrms = zeros(numel(idx_radial), 1);
        for k = 1:numel(idx_radial)
            nrms(k) = basis_norm(basis, ell, k);
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

            v(idx,:) = ang_radial'*x(mask,:);

            ind = ind + numel(idx);
            ind_ang = ind_ang+1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    v = roll_dim(v, sz_roll);
end

function v = fourier_bessel_expand(x, basis)
    basis_check_expand(x, basis);

    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    b = basis.evaluate_t(x);

    A = @(v)(basis.evaluate_t(basis.evaluate(v)));

    cg_opt.max_iter = Inf;
    cg_opt.rel_tolerance = 1e-15;
    cg_opt.verbose = false;

    v = conj_grad(A, b, cg_opt);
 
    v = roll_dim(v, sz_roll);
end

function x = fourier_bessel_expand_t(v, basis)
    basis_check_evaluate(v, basis);

    [v, sz_roll] = unroll_dim(v, 2);

    b = vol_to_vec(basis.evaluate(v));

    A = @(x)(vol_to_vec(basis.evaluate(basis.evaluate_t(vec_to_vol(x)))));

    cg_opt.max_iter = Inf;
    cg_opt.rel_tolerance = 1e-15;
    cg_opt.verbose = false;

    x = conj_grad(A, b, cg_opt);

    x = roll_dim(x, sz_roll);
end
