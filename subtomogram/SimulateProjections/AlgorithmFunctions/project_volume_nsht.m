function [recon, T, L] = project_volume_nsht(vol, angles, R, N, ns)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% spherical bessel expansion
T.start = tic;
%L = 3*R;
c = 1/2;
ss = 2;
n_r = 4*c*R; 
[coeff, basis, sample_points, L] = spherical_expansion_nsht(vol, n_r, ss);
T.expand = toc(T.start);

%% precomputation for projection
[Wp, Wm, Crm, nn, fn, r0, InnerP] = precompute_projection(ns, n_r, R, c, L,ss, basis, sample_points);
T.precompute = toc(T.start);
clear basis sample_points;

%% Rotate coefficient, restrict to xy or other planes in subtomograms
batch_size = floor(1000/ns);
recon = [];
for i = 1:ceil(N/batch_size)
    actual_batch_size = min(batch_size, N-(i-1)*batch_size);
    startidx = (i-1)*batch_size;
    CR = rotateCellFast(coeff, angles(:,startidx+1:startidx+actual_batch_size), Wp, Wm);
    %T.rotate = toc(T.start);

    % for each m, compute inner product multiply by rotated coefficients
    [coeff_pos_k] = project_coeff(CR, Crm, InnerP, L, actual_batch_size, ns, nn);
    %T.project = toc(T.start);

    % recon
    nx = 2*R+1;
    [recon2] = recon_images_FBm_2(R, nx, coeff_pos_k, 1, actual_batch_size*ns, fn, r0);%/sqrt(2*pi);
    recon = cat(3, recon, recon2);
end
T.recon = toc(T.start);
clear CR Crm InnerP;

T.end = T.recon;
T.precompute = T.precompute - T.expand;
% T.recon = T.recon - T.project;
% T.project = T.project - T.rotate;
% T.rotate = T.rotate - T.precompute;
end

