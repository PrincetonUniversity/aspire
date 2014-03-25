function [ v, v_b, kernel ,err, iter, flag] = recon3d_firm_ctf_large( projections,ctfs,defocusID,inv_rot_matrices,shifts, tol, max_it, x0 , group_size)

% Reconstruction from 2d images using FIRM with CTF correction

% Input:

%

%   projections: stack of projections, of size n x n x n_proj, where

%   n is the size of a projection, and n_proj is the number of the

%   projections.

%

%   ctfs: a stack of ctf images in Fourier space, of size n x n x n_d,

%   where n_d is the number of defocus groups.

%

%   defocusID: record the defocus group each image belongs to, an array of

%   length n_proj, the indices are in the set {1,2,...,n_d}.

%

%   shifts: the translation parameters for the projections. 

%   

%   inv_rot_matrices: a stack of inverse rotation matrices, of size 3x3xn_proj.

%

%   tol     error toleration

%   max_it  maximun number of iterations

%   x0      initial guess

%

%   

%   group_size: number of images in each group to compute the back projection v_b

%   and the kernel groupwise. Default value: 1000.

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
n_proj = size(projections, 3);

if isempty(tol), tol=1e-3; end

if isempty(max_it), max_it=100; end

if isempty(x0), x0 = zeros(n,n,n); end

if isempty(group_size), group_size=1000; end

if ~exist('group_size','var'), group_size=1000; end

v_b = 0;

kernel = 0;



for i = 1 : group_size : n_proj

 ind=i:min(i+group_size-1, n_proj);

 sub_projections=projections(:,:,ind);

 sub_defocusID=defocusID(ind);

 sub_inv_rot_matrices=inv_rot_matrices(:,:,ind);

if isempty(shifts)

    projs_fourier= FFT_projections( sub_projections);

else
    sub_shifts = shifts(ind, :);
    projs_fourier= FFT_projections( sub_projections,sub_shifts);

end

[sub_v_b,sub_kernel] = precomp_ctf( projs_fourier,sub_inv_rot_matrices, ctfs,sub_defocusID);

v_b = v_b + sub_v_b;

kernel = kernel + sub_kernel;



end

[v, err, iter, flag] = cg_recon3d(kernel, v_b, tol, max_it, x0);

end

