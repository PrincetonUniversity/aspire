function [recon, T, L] = project_volume(vol, angles, R, N, ns, input_type)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if nargin>5 && sum(input_type == 'real') > 2
    symmetric = 1;
else
    symmetric = 0;
end

%% spherical bessel expansion
T.start = tic;
%L = 3*R;
c = 1/2;
ss = 2;
n_r = 4*c*R; 
if symmetric
    [coeff, basis, sample_points, L] = spherical_expansion_ss(vol, n_r, ss);
else
    [coeff, basis, sample_points, L] = spherical_expansion(vol, n_r, ss);
end
T.expand = toc(T.start);

%% precomputation for projection

[Wp, Wm, Crm, nn, fn, r0, InnerP] = precompute_projection(ns, n_r, R, c, L,ss, basis, sample_points);

T.precompute = toc(T.start);
clear basis sample_points;

%% Rotate coefficient, restrict to xy or other planes in subtomograms
batch_size = floor(1000/ns);
recon = zeros(2*R+1, 2*R+1, batch_size*ns, ceil(N/batch_size));
angles = cat(2, angles, zeros(3, ceil(N/batch_size)*batch_size - N));
angles = reshape(angles, 3, batch_size, []);

parfor i = 1:ceil(N/batch_size)
    actual_batch_size = min(batch_size, N-(i-1)*batch_size);
    angle = angles(:,:,i);
    if actual_batch_size < batch_size
        angle = angle(:,1:actual_batch_size);
    end
    CR = rotateCellFast(coeff, angle, Wp, Wm);
    %T.rotate = toc(T.start);

    % for each m, compute inner product multiply by rotated coefficients
    [coeff_pos_k] = project_coeff(CR, Crm, InnerP, L, actual_batch_size, ns, nn);
    %T.project = toc(T.start);

    % recon
    nx = 2*R+1;
    [recon2] = recon_images_FBm_2(R, nx, coeff_pos_k, 1, actual_batch_size*ns, fn, r0);%/sqrt(2*pi);
    if actual_batch_size < batch_size
        recon2 = cat(3, recon2, zeros(2*R+1, 2*R+1, batch_size - actual_batch_size));
    end
    recon(:,:, :, i) = recon2;
end

T.recon = toc(T.start);
recon = reshape(recon, 2*R+1, 2*R+1, []);
recon = recon(:,:,1:N*ns);
clear CR Crm InnerP;

T.end = T.recon;
T.precompute = T.precompute - T.expand;
% T.recon = T.recon - T.project;
% T.project = T.project - T.rotate;
% T.rotate = T.rotate - T.precompute;
end

