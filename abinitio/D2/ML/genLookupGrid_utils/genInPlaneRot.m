
%% This function takes projection direcrions (points on the 2-sphere) and
% generates rotation matrices in SO(3). The projection direction 
% is the 3rd column and columns 1 and 2 span the perpendicular plane. 
% To properly discretize SO(3), for each projection direction we generate
% [2*pi/dtheta] "in-plane" rotations, of the plane 
% perpendicular to this direction. This is done by generating one rotation
% for each direction and then multiplying on the right by a rotation about
% the Z-axis by k*dtheta degrees, k=0...2*pi/dtheta-1. 

% Input: 
% sphere_grid = Discretized 2-sphere
% dtheta = resolution for in-plane rotations

% Output: 
% rots_grid = Set of rotations, one for each point in sphere_grid
% inplane_rotated_grid = In-plane rotations of each rotation in rots_grid
% n_inplane = number of in-plane rotations (derived from resolution dtheta) 

function [inplane_rotated_grid,rots_grid,n_inplane]=genInPlaneRot(sphere_grid,dtheta)
%% Generate one rotation for each point one the sphere
nrot=size(sphere_grid,2);
Ri2=cat(1,-(sphere_grid(2,:)),sphere_grid(1,:),zeros(1,nrot));
Ri2_norms=repmat(sqrt(sum(Ri2.^2,1)),3,1);
Ri2=Ri2./Ri2_norms;
Ri1=cross(Ri2,sphere_grid,1);
Ri1_norms=repmat(sqrt(sum(Ri1.^2,1)),3,1);
Ri1=Ri1./Ri1_norms;
rots_grid=cat(3,Ri1,Ri2,sphere_grid);
rots_grid=permute(rots_grid,[1,3,2]);

%% Generate in plane rotations. 
% Generate [2*pi/dtheta] rotations about the Z-axis. 
dtheta_vec=0:dtheta:2*pi-dtheta;
n_inplane=length(dtheta_vec);
inplane_rots=zeros(3,3,n_inplane);
cosines=cos(2*pi-dtheta_vec);
sines=sin(2*pi-dtheta_vec);
inplane_rots(1,1,:)=cosines;
inplane_rots(2,2,:)=cosines;
inplane_rots(1,2,:)=-sines;
inplane_rots(2,1,:)=sines;
inplane_rots(3,3,:)=1;

%% Generate in-plane rotations of rots_grid by multiplying on the right by 
%the in-plane rotations inplane_rots. 
inplane_rotated_grid=zeros(3,3,n_inplane,nrot);
for i=1:nrot
    inplane_rotated_grid(:,:,:,i)=...
        multiprod(rots_grid(:,:,i),inplane_rots);
end
