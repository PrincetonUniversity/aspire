function [v_b,kernel] = precomp_ctf( projs_fourier,inv_rot_matrices,ctfs,defocusID,fprecomp)
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
%   fprecomp: precomputed kernel matrix (optional). For fixed projection
%   directions.
%
% Output:
%
%   v_b: the backprojection of the Fourier slices.
%   kernel: the kernel matrix corresponding to A^*A.
%
% Lanhui Wang, Princeton University, Feb 10, 2012

precomp=0;
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

n=size(projs_fourier,1);
n_proj=size(projs_fourier,3);

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
omega=zeros(m^2*n_proj,3);
pfs=zeros(n,n,n_proj,class(projs_fourier));
weights=zeros(n,n,n_proj);
for k=1:n_proj;
    P = I * n_x(:,k)' + J * n_y(:,k)';
    omega((k-1)*m^2+1:k*m^2,:)=P;
    pf=projs_fourier(:,:,k);
    pfs(:,:,k)=pf.*ctfs(:,:,defocusID(k));
    weights(:,:,k)=ctfs(:,:,defocusID(k)).^2;
end;
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
v_b = anufft3(pfs, omega', Nd);

if precomp==0
    kernel=zeros(n*2,n*2,n*2);
    for s1=[-n,0]
        for s2=[-n,0]
            s3=0;
            idx1=s1+n+(1:n);
            idx2=s2+n+(1:n);
            idx3=s3+n+(1:n);
            shift=[-s1 -s2 -s3];
            shift = -shift+n_shift;
            kernel(idx1,idx2,idx3) = anufft3( ...
                weights.*exp(1i*(omega*shift(:))), omega', Nd);
        end
    end
    s3=-n+1;% use toepliz property
    idx3=s3+n+(1:n-1);
    cidx=-((2:2*n)-n-1)+n+1;
    kernel(2:end,2:end,idx3)=conj(kernel(cidx,cidx,-(idx3-n-1)+n+1));
    kernel=kernel(2:end,2:end,2:end);
    if exist('fprecomp','var')
        save(fprecomp,'kernel');
    end
end
end
