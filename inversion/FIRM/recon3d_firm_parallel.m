function [ v, v_b, kernel ,err, iter, flag] = recon3d_firm_parallel( projections,...
            inv_rot_matrices,shifts, tol, max_it, x0 ,fprecomp)
% Reconstruction from 2d images using FIRM
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
%
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
n=size(projections,1);
if isempty(tol), tol=1e-3; end
if isempty(max_it), max_it=100; end
if isempty(x0), x0 = zeros(n,n,n); end

if isempty(shifts)
    projs_fourier= FFT_projections( projections);
else
    projs_fourier= FFT_projections( projections,shifts);
end
if exist('fprecomp','var')
    [v_b,kernel] = precomp_parallel( projs_fourier,inv_rot_matrices,fprecomp);
else
    [v_b,kernel] = precomp_parallel( projs_fourier,inv_rot_matrices);
end
[v, err, iter, flag] = cg_recon3d(kernel, v_b, tol, max_it, x0);
end

