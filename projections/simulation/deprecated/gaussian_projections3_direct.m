function p=gaussian_projections3_direct(def,rot_matrices,n,nsteps,rmax)
%
% XXX FIX
%
% Yoel Shkolnisky, June 2007.

if nargin<4
    nsteps=500;
end

if nargin<5
    rmax=1;
end


t=linspace(-abs(rmax),abs(rmax),n);
[u,v]=meshgrid(t,t); % Do not use ndgrid to be compatible with gen_projections.
u=u(:);
v=v(:);

t0=sqrt(12)/2;
dt=2*t0/nsteps;

inv_rot_matrices = permute(rot_matrices, [2 1 3]);
p=zeros(size(u));

for k=1:numel(u)
    
    alpha = inv_rot_matrices(:,1).';
    beta  = inv_rot_matrices(:,2).';
    tau   = inv_rot_matrices(:,3).';
    
    x0=u(k)*alpha;
    y0=v(k)*beta;
    pts=repmat(x0,nsteps,1)+repmat(y0,nsteps,1)+((0:nsteps-1)*dt-t0).'*tau;
    
    p_vals=gaussian_phantom3_v2(pts,def);
    p(k)=sum(p_vals)*dt;
    
end

p=reshape(p,n,n);


