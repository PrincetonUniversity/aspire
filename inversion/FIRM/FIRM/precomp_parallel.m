function [v_b,kernel] = precomp_parallel( projs_fourier,inv_rot_matrices,fprecomp)
% Precompute the backprojection and the kernel matrix for FIRM
% Input:
%
%   projs_fourier: a stack of square Fourier slices, of size n x n x n_proj,
%   where n is the size of a square Fourier slice, and n_proj
%   is the number of the Fourier slices.
%   
%   inv_rot_matrices: a stack of inverse rotation matrices, of size 3x3xn_proj.
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
ps=matlabpool('size');
if ps==0
    matlabpool open;
end
ps=matlabpool('size');
parfor k=1:n_proj;
    P = I * n_x(:,k)' + J * n_y(:,k)';
    omega(:,k,:)=P;
end;
omega=reshape(omega,m^2*n_proj,3);
% points in the Fourier cube
ind_cub= max(omega,[],2)<fix(n/2) & min(omega,[],2)>-fix(n/2);
omega=omega(ind_cub,:);
pfs=projs_fourier(ind_cub);

%% Precompute the backprojection and the kernel matrix
% V=nufft_3d(projs_fourier,omega,precision,M);
Jd = [4 4 4];
Nd = [n n n];
Ld = [2^11 2^11 2^11];% use table for large matrix
Kd = 2 * Nd;
% 		n_shift = zeros(size(Nd)); printm 'easy: 0 shift'
n_shift = fix(Nd/2); % stress it

gam = 2*pi ./ Nd;
omega=[omega(:,1)*gam(1) omega(:,2)*gam(2) omega(:,3)*gam(3)];
ktype = 'minmax:kb';
% H = nufft_init(omega, Nd, Jd, Kd, n_shift, ktype);
v_b=zeros(n,n,n);
total_length=length(pfs);
ps=16;
core_length=ceil(total_length/ps);
parfor core=1:ps
    core_ind= ((core-1)*core_length+1):min(core*core_length,total_length);
    core_ind=core_ind(:);
    st = nufft_init(omega(core_ind,:), Nd, Jd, Kd, n_shift, 'table', Ld, ktype);
    v_b_temp = nufft_adj(pfs(core_ind), st);
    v_b=v_b+v_b_temp;
    kernel_temp=zeros(n*2,n*2,n*2);
    if precomp==0
        for s1=[-n,0]
            for s2=[-n,0]
                s3=0;
                idx1=s1+n+(1:n);
                idx2=s2+n+(1:n);
                idx3=s3+n+(1:n);
                st1=st;
                shift=[-s1 -s2 -s3];
                st1.phase_shift = exp(1i * (st.om * shift(:)));
                %         kernel(idx1,idx2,idx3) = nufft_adj(ones(size(ind_cub)), st1);
                kernel_temp(idx1,idx2,idx3) = nufft_adj(ones(size(core_ind)), st1);
            end
        end
        kernel=kernel+kernel_temp;
    end
end
matlabpool close;
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