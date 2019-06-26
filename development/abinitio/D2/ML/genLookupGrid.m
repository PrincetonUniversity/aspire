
%Input: N is size of grid, best choose N=k^2
function [debug_data]=genLookupGrid(N,eq_filter_angle,dtheta,s)
%% Generate uniform grid on sphere with Saaf-Kuijlaars

if nargin>3
    rng(s);
end

sphere_grid=SaffKuijlaars(N);
octant1_idx=(sphere_grid(:,1)>0).*(sphere_grid(:,2)>0).*(sphere_grid(:,3)>0);
octant2_idx=(sphere_grid(:,1)<0).*(sphere_grid(:,2)>0).*(sphere_grid(:,3)>0);
sphere_grid2=sphere_grid(octant2_idx==1,:);
sphere_grid=sphere_grid(octant1_idx==1,:);

%DEBUG
% sphere_grid2=sphere_grid;
% sphere_grid2(:,1:2)=sphere_grid2(:,[2,1]);
% sphere_grid2(:,1)=-sphere_grid2(:,1);

nrot=size(sphere_grid,1);
nrot2=size(sphere_grid2,1);

%% Mark points on equators
%project each vector onto xy,xz,yz planes and measure angular distance
sphere_grid=sphere_grid'; 
sphere_grid2=sphere_grid2';
eq_idx=markEquators(sphere_grid,eq_filter_angle);
eq_idx2=markEquators(sphere_grid2,eq_filter_angle);

%% Determine in-plane angular resolution to be ~ sphere grid resolution
if nargin<3
    pi_devisors=[1,2,3,4,5,6,9,10,12,15,20,30,60,90,180];    
    mean_angular_dist=getGridMeanDist(sphere_grid);
    [~,I]=min(abs(pi_devisors-mean_angular_dist));
    dtheta=pi_devisors(I);
else
    mean_angular_dist=dtheta;
end
ntheta=360/dtheta;
dtheta=dtheta*pi/180;
%% Generate random roatations for testing by chooosing Ri(:,2) randomly
[inplane_rotated_grid,rots_grid,n_inplane]=genInPlaneRot(sphere_grid,dtheta);
[inplane_rotated_grid2,rots_grid2,n_inplane2]=genInPlaneRot(sphere_grid2,dtheta);

%% Generate all relative rotations in octant 1
inplane_rotated_grid=single(inplane_rotated_grid);
eq_idx=single(logical(eq_idx));
eq2eq_Rij_table=logical(triu(ones(nrot,nrot)-eq_idx*eq_idx',1));
[unique_pairs,~]=removeD2Ambiguities(sphere_grid);
unique_pairs=logical(unique_pairs.*eq2eq_Rij_table);
npairs=sum(unique_pairs(:));
h=waitbar(0,'Generating relative rotations octant 1...');
idx=0;
idx_vec=1:nrot;

cls=zeros(2,4,0.5*ntheta,ntheta,2*npairs);

for i=1:nrot-1
   
    unique_pairs_i=idx_vec(unique_pairs(i,:));
    n2=sum(unique_pairs_i);
    if n2==0
        continue;
    end
    for j=unique_pairs_i      
        idx=idx+1;
        Rjs=inplane_rotated_grid(:,:,1:0.5*n_inplane,j);
        Ris=permute(inplane_rotated_grid(:,:,:,i),[1,2,4,3]);
        Rijs=multiprod(multi_transpose(Rjs),Ris);
        
        cls(1,1,:,:,npairs+idx)=atan2(Rijs(1,3,:,:),-Rijs(2,3,:,:));
        cls(2,1,:,:,npairs+idx)=atan2(-Rijs(3,1,:,:),Rijs(3,2,:,:));
            cls(1,1,:,:,idx)=atan2(Rijs(3,1,:,:),-Rijs(3,2,:,:));
            cls(2,1,:,:,idx)=atan2(-Rijs(1,3,:,:),Rijs(2,3,:,:));        
        
        gRjs=Rjs;
        gRjs(2:3,:,:)=-gRjs(2:3,:,:);
        Rijs=multiprod(multi_transpose(gRjs),Ris);
        
        cls(1,2,:,:,npairs+idx)=atan2(Rijs(1,3,:,:),-Rijs(2,3,:,:));
        cls(2,2,:,:,npairs+idx)=atan2(-Rijs(3,1,:,:),Rijs(3,2,:,:));
        cls(1,2,:,:,idx)=atan2(Rijs(3,1,:,:),-Rijs(3,2,:,:));
        cls(2,2,:,:,idx)=atan2(-Rijs(1,3,:,:),Rijs(2,3,:,:)); 
        
        gRjs=Rjs;
        gRjs([1,3],:,:)=-gRjs([1,3],:,:);
        Rijs=multiprod(multi_transpose(gRjs),Ris);
        
        cls(1,3,:,:,npairs+idx)=atan2(Rijs(1,3,:,:),-Rijs(2,3,:,:));
        cls(2,3,:,:,npairs+idx)=atan2(-Rijs(3,1,:,:),Rijs(3,2,:,:));
        cls(1,3,:,:,idx)=atan2(Rijs(3,1,:,:),-Rijs(3,2,:,:));
        cls(2,3,:,:,idx)=atan2(-Rijs(1,3,:,:),Rijs(2,3,:,:)); 
        
        gRjs=Rjs;
        gRjs(1:2,:,:)=-gRjs(1:2,:,:);
        Rijs=multiprod(multi_transpose(gRjs),Ris);

        cls(1,4,:,:,npairs+idx)=atan2(Rijs(1,3,:,:),-Rijs(2,3,:,:));
        cls(2,4,:,:,npairs+idx)=atan2(-Rijs(3,1,:,:),Rijs(3,2,:,:));
        cls(1,4,:,:,idx)=atan2(Rijs(3,1,:,:),-Rijs(3,2,:,:));
        cls(2,4,:,:,idx)=atan2(-Rijs(1,3,:,:),Rijs(2,3,:,:)); 
    end
    waitbar(idx/npairs);
end
close(h);

%% Generate all relative rotations in octant 2
inplane_rotated_grid2=single(inplane_rotated_grid2);
eq_idx2=single(logical(eq_idx2));
eq2eq_Rij_table=logical(ones(nrot,nrot2)-eq_idx*eq_idx2');
%[unique_pairs,~]=removeD2Ambiguities(sphere_grid2);
unique_pairs12=eq2eq_Rij_table;
npairs12=sum(unique_pairs12(:));
h=waitbar(0,'Generating relative rotations, octant 1 to octant 2...');
Rijs_grid=[];
idx=0;
idx_vec=1:nrot2;

cls2=zeros(2,4,0.5*ntheta,ntheta,2*npairs12);

for i=1:nrot
   
    unique_pairs_i=idx_vec(unique_pairs12(i,:));
    n2=sum(unique_pairs_i);
    if n2==0
        continue;
    end
    for j=unique_pairs_i      
        idx=idx+1;
        Rjs=inplane_rotated_grid2(:,:,1:0.5*n_inplane2,j);
        Ris=permute(inplane_rotated_grid(:,:,:,i),[1,2,4,3]);
        Rijs=multiprod(multi_transpose(Rjs),Ris);
        
        cls2(1,1,:,:,npairs12+idx)=atan2(Rijs(1,3,:,:),-Rijs(2,3,:,:));
        cls2(2,1,:,:,npairs12+idx)=atan2(-Rijs(3,1,:,:),Rijs(3,2,:,:));
        cls2(1,1,:,:,idx)=atan2(Rijs(3,1,:,:),-Rijs(3,2,:,:));
        cls2(2,1,:,:,idx)=atan2(-Rijs(1,3,:,:),Rijs(2,3,:,:));        
        
        gRjs=Rjs;
        gRjs(2:3,:,:)=-gRjs(2:3,:,:);
        Rijs=multiprod(multi_transpose(gRjs),Ris);
        
        cls2(1,2,:,:,npairs12+idx)=atan2(Rijs(1,3,:,:),-Rijs(2,3,:,:));
        cls2(2,2,:,:,npairs12+idx)=atan2(-Rijs(3,1,:,:),Rijs(3,2,:,:));
        cls2(1,2,:,:,idx)=atan2(Rijs(3,1,:,:),-Rijs(3,2,:,:));
        cls2(2,2,:,:,idx)=atan2(-Rijs(1,3,:,:),Rijs(2,3,:,:)); 
        
        gRjs=Rjs;
        gRjs([1,3],:,:)=-gRjs([1,3],:,:);
        Rijs=multiprod(multi_transpose(gRjs),Ris);
        
        cls2(1,3,:,:,npairs12+idx)=atan2(Rijs(1,3,:,:),-Rijs(2,3,:,:));
        cls2(2,3,:,:,npairs12+idx)=atan2(-Rijs(3,1,:,:),Rijs(3,2,:,:));
        cls2(1,3,:,:,idx)=atan2(Rijs(3,1,:,:),-Rijs(3,2,:,:));
        cls2(2,3,:,:,idx)=atan2(-Rijs(1,3,:,:),Rijs(2,3,:,:)); 
        
        gRjs=Rjs;
        gRjs(1:2,:,:)=-gRjs(1:2,:,:);
        Rijs=multiprod(multi_transpose(gRjs),Ris);

        cls2(1,4,:,:,npairs12+idx)=atan2(Rijs(1,3,:,:),-Rijs(2,3,:,:));
        cls2(2,4,:,:,npairs12+idx)=atan2(-Rijs(3,1,:,:),Rijs(3,2,:,:));
        cls2(1,4,:,:,idx)=atan2(Rijs(3,1,:,:),-Rijs(3,2,:,:));
        cls2(2,4,:,:,idx)=atan2(-Rijs(1,3,:,:),Rijs(2,3,:,:)); 
    end
    waitbar(idx/npairs12);
end

close(h);
clearvars Ris Rjs

%% Generate fourier ray angles (indexes) from Rijs_grid
disp('Generating common lines...');
cls=mod(cls+2*pi,2*pi); %make all angles non negative 
cls=cls*180/pi; %convert to degrees 
cls=reshape(cls,2,4*0.5*ntheta^2*npairs*2)';
[cls_lookup]=getClsFromRijs3(cls);
%cls_lookup=reshape(cls_lookup,4,0.5*ntheta*ntheta*npairs*2);
cls=uint16(cls);
cls11=reshape(cls(:,1),4,1,0.5*ntheta^2*npairs*2);
cls12=reshape(cls(:,2),4,1,0.5*ntheta^2*npairs*2);
cls=cat(2,cls11,cls12);

%% Generate fourier ray angles (indexes) from Rijs_grid octant 2
disp('Generating common lines...');
cls2=mod(cls2+2*pi,2*pi); %make all angles non negative 
cls2=cls2*180/pi; %convert to degrees 
cls2=reshape(cls2,2,4*0.5*ntheta^2*npairs12*2)';
[cls2_lookup]=getClsFromRijs3(cls2);
%cls_lookup=reshape(cls_lookup,4,0.5*ntheta*ntheta*npairs12*2);
cls2=uint16(cls2);
cls11=reshape(cls2(:,1),4,1,0.5*ntheta^2*npairs12*2);
cls12=reshape(cls2(:,2),4,1,0.5*ntheta^2*npairs12*2);
cls2=cat(2,cls11,cls12);
%cls_lookup=[cls_lookup,cls2_lookup];
%% Generate rots_grid statistics for debuging
% latitude line angles and longitude values with respect to all symmetry 
% axes, angles between projection directions and 2 angles of spherical
% coordinates representation. 
[grid_stats,pol_rep]=generateGridStats(rots_grid);

%% Finalize
debug_data=struct('rots_grid',rots_grid,'cls_lookup',cls_lookup,...
    'Rijs_grid',Rijs_grid,'grid_stats',grid_stats,'nrot',nrot,...
    'ntheta',ntheta,'dtheta',dtheta,'npairs',npairs,'eq_filter_angle',...
    eq_filter_angle,'mean_angular_dist',mean_angular_dist,'pol_rep',...
    pol_rep,'cls',cls,'inplane_rotated_grid',inplane_rotated_grid,...
    'unique_pairs',unique_pairs,'cls2_lookup',cls2_lookup,'cls2',cls2,...
    'nrot2',nrot2,'unique_pairs12',unique_pairs12,'inplane_rotated_grid2',...
    inplane_rotated_grid2,'rots_grid2',rots_grid2,'npairs12',npairs12,...
    'eq_idx',eq_idx,'eq_idx2',eq_idx2);
end