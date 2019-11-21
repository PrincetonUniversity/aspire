function [ projs ] = Gaussian_gen_projections( R,K,refq, n_theta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu1 = [1 0 0]; 
mu2 = [0 1 0]; 
mu3 = [0 0 1];
Sigma1 = [1 2 2]*16; 
Sigma2 = [2 3 2]*16; 
Sigma3 = [3 3 4]*16;

% sample volume in real domain and NUFFT transform to Fourier domain
nx = 2*R+1;
x1 = -R:1:R;
[X,Y,Z] = meshgrid(x1,x1,x1);

F1 = mvnpdf([X(:) Y(:) Z(:)],mu1,Sigma1);
F2 = mvnpdf([X(:) Y(:) Z(:)],mu2,Sigma2);
F3 = mvnpdf([X(:) Y(:) Z(:)],mu3,Sigma3);
F3d = reshape(F1+F1+F3,nx,nx,nx); 

% samples in Fourier space
c = pi/2;
dtheta = 2*pi/n_theta;
dr = c/n_theta;
thetas = 0:dtheta:2*pi-dtheta;
rs = 0:dr:c-dr;
x = rs'*cos(thetas);
y = rs'*sin(thetas);
z = zeros(n_theta);
omega = [x(:) y(:) z(:)]';
projs = zeros(n_theta,n_theta,K);

parfor i = 1:K;
    rot = q_to_rot(refq(:,i));
    omega_rot = rot'*omega;
    projection = nufft3d2(n_theta^2,omega_rot(1,:),omega_rot(2,:),...
        omega_rot(3,:),-1,1e-15,nx,nx,nx,F3d);
    projs(:,:,i) = reshape(projection,n_theta,n_theta);
end

end

