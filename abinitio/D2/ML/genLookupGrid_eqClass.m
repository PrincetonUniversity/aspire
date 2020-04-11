
%Input: N is size of grid, best choose N=k^2
function [lookup_data]=genLookupGrid_eqClass(N,eq_filter_angle,dtheta,s)

% Use pre-defined seed to replicate results. 
if nargin>3
    rng(s);
end

%% Generate uniform grid on sphere with Saaf-Kuijlaars and take one quarter
 % of sphere because of D2 symmetry redundancy. sphere_grid is the
 % upper part of the quarter (above XY plane and sphere_grid2 is the lower
 % part (below XY plane). 
sphere_grid=SaffKuijlaars(N);
octant1_idx=(sphere_grid(:,1)>0).*(sphere_grid(:,2)>0).*(sphere_grid(:,3)>0);
octant2_idx=(sphere_grid(:,1)>0).*(sphere_grid(:,2)>0).*(sphere_grid(:,3)<0);
sphere_grid2=sphere_grid(octant2_idx==1,:);
sphere_grid=sphere_grid(octant1_idx==1,:);

%% Mark Equator directions
%  Common lines between projection directions which are perpendicular to 
%  symmetry axes (equator images) have common line degeneracies. Two images
%  taken from directions on the same great circle which is perpendicular to
%  some symmetry axis only have 2 common lines instead of 4, and must be 
%  treated separately. 
%  We detect such directions by taking a strip of radius
%  eq_filter_angle about the 3 great circles perpendicular to the symmetry
%  axes of D2 (i.e to X,Y and Z axes). 
sphere_grid=sphere_grid'; 
sphere_grid2=sphere_grid2';
[eq_idx,eq_class]=markEquators(sphere_grid,eq_filter_angle);
[eq_idx2,eq_class2]=markEquators(sphere_grid2,eq_filter_angle);

%% Mark Top view directions
%  A Top view projection image is taken from the direction of one of the
%  symmetry axes. Since all symmetry axes of D2 molecules are perpendicular
%  this means that such an image is an equator with repect to both symmetry
%  axes which are perpendicular to the direction of the symmetry axis from
%  which the image was made, e.g. if the image was formed by projecting in
%  the direction of the X (symmetry) axis, then it is an equator with
%  respect to both Y and Z symmetry axes (it's direction is the
%  interesection of 2 great circles perpendicular to Y and Z axes). 
%  Such images have severe degeneracies. A pair of Top View images (taken
%  from different directions or a Top View and equator image only have a 
%  single common line. A top view and a regular non-equator image only have
%  two common lines. 
non_tv1=eq_class<=3;
non_tv2=eq_class2<=3;
sphere_grid=sphere_grid(:,non_tv1==1);
sphere_grid2=sphere_grid2(:,non_tv2==1);
eq_idx=eq_idx(non_tv1==1);
eq_idx2=eq_idx2(non_tv2==1);
eq_class=eq_class(non_tv1==1);
eq_class2=eq_class2(non_tv2==1);
nrot=size(sphere_grid,2);
nrot2=size(sphere_grid2,2);

% %DEBUG
% eq_idx=zeros(size(eq_idx));
% eq_idx2=zeros(size(eq_idx2));
% eq_class=zeros(size(eq_class));
% eq_class2=zeros(size(eq_class2));

%% Determine the resolution for in-plane resolution
%  If not given by user then in-plane resolution ~ sphere grid resolution
%  we compute mean angular distance between points on the discretization of
%  the sphere and round it. 
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

%% Generate all relative rotations candtates for ML in octant 1 (above XY plane)
%  Do everything in single precision to save space. 
inplane_rotated_grid=single(inplane_rotated_grid);
eq_idx=single(eq_idx);

%  Generate upper triangular table of indicators of all pairs which are not
%  equators with respect to the same symmetry axis (named unique_pairs). 
eq_table_idx=eq_idx*eq_idx';
eq2eq_Rij_table=triu(ones(nrot,nrot)-...
    eq_table_idx.*(abs(round(eq_class-eq_class'))==0),1);
%eq2eq_Rij_table11=eq2eq_Rij_table;
%eq2eq_Rij_table11=logical(eq2eq_Rij_table);
eq2eq_Rij_table11=logical(eq2eq_Rij_table);

npairs=sum(eq2eq_Rij_table11(:));
%h=waitbar(0,'Generating relative rotations octant 1...');
idx=0;
idx_vec=1:nrot;
cls=zeros(2,4,0.5*ntheta,ntheta,2*npairs);
for i=1:nrot-1
   
    unique_pairs_i=idx_vec(eq2eq_Rij_table11(i,:));
    n2=sum(unique_pairs_i);
    if n2==0
        continue;
    end
    for j=unique_pairs_i 
        
        % Compute relative rotations candidates
        idx=idx+1;
        Rjs=inplane_rotated_grid(:,:,1:0.5*n_inplane,j);
        Ris=permute(inplane_rotated_grid(:,:,:,i),[1,2,4,3]);
        Rijs=multiprod(multi_transpose(Rjs),Ris);
        
        % Compute induced common lines 
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
    %waitbar(idx/npairs);
end
%close(h);

%% Generate all relative rotations in octant 2
%  This is the same as with octant 1 but relative rotations are computed
%  for directions where one is taken in octant 1 and the other in octant 2
inplane_rotated_grid2=single(inplane_rotated_grid2);
eq_idx2=single(eq_idx2);

eq_table_idx=eq_idx*eq_idx2';%places of eq pairs in the table
eq2eq_Rij_table=ones(nrot,nrot2)-...
    eq_table_idx.*(abs(round(eq_class-eq_class2'))==0);
eq2eq_Rij_table12=logical(eq2eq_Rij_table);
%unique_pairs12=logical(eq2eq_Rij_table);

npairs12=sum(eq2eq_Rij_table12(:));
%h=waitbar(0,'Generating relative rotations, octant 1 to octant 2...');
Rijs_grid=[];
idx=0;
idx_vec=1:nrot2;
cls2=zeros(2,4,0.5*ntheta,ntheta,2*npairs12);

for i=1:nrot
   
    unique_pairs_i=idx_vec(eq2eq_Rij_table12(i,:));
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
    %waitbar(idx/npairs12);
end

%close(h);
clearvars Ris Rjs

log_message('Generating lookup data for common lines estimation');

%% Generate fourier ray angles indices into correlations tables
%  for common lines between projection directions in octant 1 
cls=mod(cls+2*pi,2*pi); %make all angles non negative 
cls=cls*180/pi; %convert to degrees 
cls=reshape(cls,2,4*0.5*ntheta^2*npairs*2)';
[cls_lookup]=getClsFromRijs3(cls);
cls=uint16(cls); % cls are now indices between 1 and 180*360, save space
cls11=reshape(cls(:,1),4,1,0.5*ntheta^2*npairs*2);
cls12=reshape(cls(:,2),4,1,0.5*ntheta^2*npairs*2);
cls=cat(2,cls11,cls12);

%% Generate fourier ray angles indices into correlations tables
%  for common lines between projection directions in octant 1 and octant 2
%disp('Generating common lines...');
cls2=mod(cls2+2*pi,2*pi); %make all angles non negative 
cls2=cls2*180/pi; %convert to degrees 
cls2=reshape(cls2,2,4*0.5*ntheta^2*npairs12*2)';
[cls2_lookup]=getClsFromRijs3(cls2);
cls2=uint16(cls2);
cls11=reshape(cls2(:,1),4,1,0.5*ntheta^2*npairs12*2);
cls12=reshape(cls2(:,2),4,1,0.5*ntheta^2*npairs12*2);
cls2=cat(2,cls11,cls12);

%% Generate rots_grid statistics for debuging
% latitude line angles and longitude values with respect to all symmetry 
% axes, angles between projection directions and 2 angles of spherical
% coordinates representation. 
% [grid_stats,pol_rep]=generateGridStats(rots_grid);

%% Register lookup data
lookup_data=struct('rots_grid',rots_grid,'cls_lookup',cls_lookup,...
    'Rijs_grid',Rijs_grid,'nrot',nrot,...
    'ntheta',ntheta,'dtheta',dtheta,'npairs',npairs,'eq_filter_angle',...
    eq_filter_angle,'mean_angular_dist',mean_angular_dist,'cls',cls,'inplane_rotated_grid',inplane_rotated_grid,...
    'cls2_lookup',cls2_lookup,'cls2',cls2,...
    'nrot2',nrot2,'inplane_rotated_grid2',...
    inplane_rotated_grid2,'rots_grid2',rots_grid2,'npairs12',npairs12,...
    'eq_idx',eq_idx,'eq_idx2',eq_idx2,'eq2eq_Rij_table11',eq2eq_Rij_table11,...
    'eq2eq_Rij_table12',eq2eq_Rij_table12);

%% DELETE AFTER VERIFICATION
% lookup_data=struct('rots_grid',rots_grid,'cls_lookup',cls_lookup,...
%     'Rijs_grid',Rijs_grid,'grid_stats',grid_stats,'nrot',nrot,...
%     'ntheta',ntheta,'dtheta',dtheta,'npairs',npairs,'eq_filter_angle',...
%     eq_filter_angle,'mean_angular_dist',mean_angular_dist,'pol_rep',...
%     pol_rep,'cls',cls,'inplane_rotated_grid',inplane_rotated_grid,...
%     'unique_pairs',unique_pairs,'cls2_lookup',cls2_lookup,'cls2',cls2,...
%     'nrot2',nrot2,'unique_pairs12',unique_pairs12,'inplane_rotated_grid2',...
%     inplane_rotated_grid2,'rots_grid2',rots_grid2,'npairs12',npairs12,...
%     'eq_idx',eq_idx,'eq_idx2',eq_idx2,'eq2eq_Rij_table11',logical(eq2eq_Rij_table11),...
%     'eq2eq_Rij_table12',logical(eq2eq_Rij_table12));
end