function [ v, v_b, kernel ,err, iter, flag] = recon3d_firm_outofcore...
    ( projs_fname, inv_rot_matrices,shifts, tol, max_it, x0, bufsize)
% XXX FIX
%
% RECON3D_FIRM_OUTOFCORE    3D Reconstruction from 2D images using FIRM
%
% Input:
%
%   projections: stack of projections, of size n x n x n_proj, where
%   n is the size of a projection, and n_proj is the number of the
%   projections.
%
%   shifts: the translation parameters for the projections.
%
%   inv_rot_matrices: a stack of inverse rotation matrices, of size 3x3xn_proj.
%
%   tol     error toleration
%   max_it  maximun number of iterations
%   x0      initial guess
%   bufsize Size of projections buffer. Determine the number of
%            projections to hold in memory.
% Output:
%
%   v: the reconstructed volume, of size nxnxn.
%   v_b: the backprojection of the Fourier slices.
%   kernel: the kernel matrix corresponding to A^*A.
%
%   err     estimated error
%   iter    number of iterations
%   flag    INTEGER: 0 = solution found to tolerance
%                    1 = no convergence given max_it
%
% Lanhui Wang, Princeton University, Feb 10, 2012

% instack=imagestackReader(instackname,100);
% n=instack.dim(3);

projreader=imagestackReader(projs_fname);
n_proj=projreader.dim(3);
n=projreader.dim(1);

if isempty(tol), tol=1e-3; end
if isempty(max_it), max_it=100; end
if isempty(x0), x0 = zeros(n,n,n); end

if ~exist('bufsize','var')
    group_size=1000;
else
    BYTES_PER_FLOAT=4; % We are using single precision.
    group_size=max(ceil(bufsize/(BYTES_PER_FLOAT*n*n)),10); % User at least 100 images.
end

v_b = 0;
kernel = 0;

for i = 1 : group_size : n_proj
    disp(i);
    ind=i:min(i+group_size-1, n_proj);
    sub_projections=projreader.getImage(ind);
    sub_inv_rot_matrices=inv_rot_matrices(:,:,ind);
    
    if isempty(shifts)
        projs_fourier= FFT_projections( sub_projections);
    else
        sub_shifts = shifts(ind, :);
        projs_fourier= FFT_projections( sub_projections,sub_shifts);
    end
    
    [sub_v_b,sub_kernel] = precomp( projs_fourier,sub_inv_rot_matrices);
    v_b = v_b + sub_v_b;
    kernel = kernel + sub_kernel;
    
end

[v, err, iter, flag] = cg_recon3d(kernel, v_b, tol, max_it, x0);

end

