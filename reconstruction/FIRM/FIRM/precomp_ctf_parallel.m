function [v_b,kernel] = precomp_ctf_parallel( projs_fourier,...
    inv_rot_matrices,ctfs,defocusID,fprecomp)
% Precompute the backprojection and the kernel matrix for FIRM
% Input:
%
%   projs_fourier: a stack of square Fourier slices, of size n x n x n_proj,
%   where n is the size of a square Fourier slice, and n_proj
%   is the number of the Fourier slices.
%   
%   inv_rot_matrices: a stack of inverse rotation matrices, of size 3x3xn_proj.
%
%   ctfs: a stack of ctf images in Fourier space, of size n x n x n_d,
%   where n_d is the number of defocus groups.
%
%   defocusID: record the defocus group each image belongs to, an array of
%   length n_proj, the indices are in the set {1,2,...,n_d}.
%
% Output:
%
%   v_b: the backprojection of the Fourier slices.
%   kernel: the kernel matrix corresponding to A^*A.
%
% Lanhui Wang, Princeton University, Feb 10, 2012
n=size(projs_fourier,1);
n_proj=size(projs_fourier,3);

precomp=0;
kernel=zeros(n*2,n*2,n*2);
if exist('fprecomp','var')
    [pathstr, name, ext] = fileparts(fprecomp);
    if isempty(ext)
        fprecomp=[fprecomp,'.mat'];
    end
    if exist(fprecomp,'file')
        fprintf('Loading precomputed kernel...');
        load(fprecomp); % loads projections, noisy_projections, shifts, q,and N.
        precomp=1;
        fprintf('Finshed!\n');
    end
end
%% Compute the 3D locations of the Fourier slices
%%%%%% Calculate great circles over S^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The circle is uniquely described by the orthogonal system (n_x,n_y,n_z)
% n_z is perpendicular to the plane of the circle,
% which is given by the two vectors n_x and n_y
%
% n_x, n_y, n_z are the image of the unit vectors in the x, y, z
% directions under the inverse rotation

n_x(:,:) = inv_rot_matrices(:,1,:); 
n_y(:,:) = inv_rot_matrices(:,2,:);
%n_z(:,:) = inv_rot_matrices(:,3,:);  % not used - just for completeness

if mod(n,2)==0
    range=-fix(n/2):fix(n/2)-1;
else
    range=-fix(n/2):fix(n/2);
end

m=length(range);
[I,J]=meshgrid(range,range);
I=I(:);
J=J(:);
omega=zeros(m^2,n_proj,3);
ctfs=ctfs(:,:,defocusID);
pfs=zeros(n,n,n_proj);
weights=zeros(n,n,n_proj);

parfor k=1:n_proj;
    P = I * n_x(:,k)' + J * n_y(:,k)';
    omega(:,k,:)=P;
    pf=projs_fourier(:,:,k);
    pfs(:,:,k)=pf.*ctfs(:,:,k);
    weights(:,:,k)=ctfs(:,:,k).^2;
end;
omega=reshape(omega,m^2*n_proj,3);
% points in the Fourier cube
ind_cub= max(omega,[],2)<fix(n/2) & min(omega,[],2)>-fix(n/2);
omega=omega(ind_cub,:);
pfs=pfs(ind_cub);
weights=weights(ind_cub);
%% Precompute the backprojection and the kernel matrix
Nd = [n n n];
n_shift = fix(Nd/2); % stress it

gam = 2*pi ./ Nd;
omega=[omega(:,1)*gam(1) omega(:,2)*gam(2) omega(:,3)*gam(3)];
v_b=zeros(n,n,n);
total_length=length(pfs);
core_length=round(total_length/ps);
parfor core=1:ps
    core_ind= ((core-1)*core_length+1):min(core*core_length,total_length);
    core_ind=core_ind(:);
    v_b_temp = anufft3(pfs(core_ind), omega(core_ind,:)', Nd);
    v_b=v_b+v_b_temp;
    if precomp==0
        kernel_temp=zeros(n*2,n*2,n*2);
        for s1=[-n,0]
            for s2=[-n,0]
                s3=0;
                idx1=s1+n+(1:n);
                idx2=s2+n+(1:n);
                idx3=s3+n+(1:n);
                shift=[-s1 -s2 -s3];
                shift = -shift+n_shift;
                kernel_temp(idx1,idx2,idx3) = anufft3( ...
                   weights(core_ind).*exp(1i*(omega(core_ind,:)*shift(:))), ...
                   omega(core_ind,:)', Nd);
            end
        end
        kernel=kernel+kernel_temp;
    end
end

if precomp==0
    s3=-n+1;% use toepliz property
    idx3=s3+n+(1:n-1);
    cidx=-((2:2*n)-n-1)+n+1;
    kernel(2:end,2:end,idx3)=conj(kernel(cidx,cidx,-(idx3-n-1)+n+1));
    kernel=kernel(2:end,2:end,2:end);
end
if exist('fprecomp','var')
    save(fprecomp,'kernel');
end
end
