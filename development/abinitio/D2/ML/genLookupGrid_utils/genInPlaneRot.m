
function [inplane_rotated_grid,rots_grid,n_inplane]=genInPlaneRot(sphere_grid,dtheta,s)
%% Generate random roatations for testing by chooosing Ri(:,2) randomly
nrot=size(sphere_grid,2);
Ri2=cat(1,-(sphere_grid(2,:)),sphere_grid(1,:),zeros(1,nrot));
Ri2_norms=repmat(sqrt(sum(Ri2.^2,1)),3,1);
Ri2=Ri2./Ri2_norms;
Ri1=cross(Ri2,sphere_grid,1);
Ri1_norms=repmat(sqrt(sum(Ri1.^2,1)),3,1);
Ri1=Ri1./Ri1_norms;
rots_grid=cat(3,Ri1,Ri2,sphere_grid);
rots_grid=permute(rots_grid,[1,3,2]);

if exist('s','var')
    rng(s);
end
rand_ang=rand(nrot,1)*2*pi;
inplane_rots=zeros(3,3,nrot);
cosines=cos(2*pi-rand_ang);
sines=sin(2*pi-rand_ang);
inplane_rots(1,1,:)=cosines;
inplane_rots(2,2,:)=cosines;
inplane_rots(1,2,:)=-sines;
inplane_rots(2,1,:)=sines;
inplane_rots(3,3,:)=1;
%rots_grid=multiprod(rots_grid,inplane_rots);

%% Rotate grid to create a lookup grid
%Generate rotation matrices
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

%Generate all in-plane rotations of rotations grid
inplane_rotated_grid=zeros(3,3,n_inplane,nrot);
for i=1:nrot
    inplane_rotated_grid(:,:,:,i)=...
        multiprod(rots_grid(:,:,i),inplane_rots);
end
