
function [scls_lookup_data]=genSelfClsDn2(lookup_data,eqf)

if nargin>1
    eq_filter_angles=eqf;
else
    eq_filter_angles=lookup_data.eq_filter_angles;
end

n=lookup_data.n;
L=360;
rots_grid1=lookup_data.rots_grid;
rots_grid2=lookup_data.rots_grid2;
sphere_grid1=squeeze(rots_grid1(:,3,:));
sphere_grid2=squeeze(rots_grid2(:,3,:));
[eq_idx,eq_class1]=markEquatorsDn2(sphere_grid1,eq_filter_angles,n);
[eq_idx2,eq_class2]=markEquatorsDn2(sphere_grid2,eq_filter_angles,n);
rots_grid1=lookup_data.inplane_rotated_grid;
rots_grid2=lookup_data.inplane_rotated_grid2;
dtheta=lookup_data.dtheta;
ntheta=lookup_data.ntheta;
[g1,g2]=gns(n);
gs=cat(3,g1(:,:,2:end),g2); %First rotation is identity

s=size(rots_grid1);
Riis1=zeros([s(1:2),2*n-1,s(3:4)]);
scls1=zeros([2,2*n-1,s(3:4)]);
nrot=s(4);
for i=1:nrot
    Ris=rots_grid1(:,:,:,i);
    for j=1:(2*n-1)
        gRis=multiprod(gs(:,:,j),Ris);
        %gRis(1:2,:,:)=-gRis(1:2,:,:);
        Riis1(:,:,j,:,i)=multiprod(permute(Ris,[2,1,3]),gRis);
        scls1(1,j,:,i)=atan2(Riis1(3,1,j,:,i),-Riis1(3,2,j,:,i));
        scls1(2,j,:,i)=atan2(-Riis1(1,3,j,:,i),Riis1(2,3,j,:,i));
    end   
end

s=size(rots_grid2);
Riis2=zeros([s(1:2),2*n-1,s(3:4)]);
scls2=zeros([2,2*n-1,s(3:4)]);
nrot2=s(4);
for i=1:nrot2
    Ris=rots_grid2(:,:,:,i);
    for j=1:2*n-1
        gRis=multiprod(gs(:,:,j),Ris);
        %gRis(1:2,:,:)=-gRis(1:2,:,:);
        Riis2(:,:,j,:,i)=multiprod(multi_transpose(Ris),gRis);
        scls2(1,j,:,i)=atan2(Riis2(3,1,j,:,i),-Riis2(3,2,j,:,i));
        scls2(2,j,:,i)=atan2(-Riis2(1,3,j,:,i),Riis2(2,3,j,:,i));
    end
end

%% Prepare common line coordinates
scls1=mod(scls1,2*pi); %make all angles non negative 
scls2=mod(scls2,2*pi); %make all angles non negative
scls1=round(scls1*180/pi);
scls2=round(scls2*180/pi);


%% Create candidate common line linear indices lists for equators
nonTv_eq_idx1=find(eq_class1>0 & eq_class1<4)';
n_eq1=length(nonTv_eq_idx1);
nonTv_eq_idx2=find(eq_class2>0 & eq_class2<4)';
n_eq2=length(nonTv_eq_idx2);
count_eq=0;

%Take care of degeneracies, for degenerate lines put the common line
%coordinates to be (1,1) so that the correlations are 1 and don't effect
%the score. 
if mod(n,2)==0 
    s1=size(rots_grid1);
    s2=size(rots_grid2);
    n_eq_idx1=sum(eq_idx==1);
    n_eq_idx2=sum(eq_idx2==1);
    scls1(:,n+1,:,[eq_class1;eq_class2]==1)=repmat([1,1],s1(3),n_eq_idx1);
    scls2(:,n+1,:,[eq_class1;eq_class2]==1)=repmat([1,1],s2(3),n_eq_idx2);
    n_eq_idx1=sum(eq_idx==2);
    n_eq_idx2=sum(eq_idx2==2);
    scls1(:,n+1,:,[eq_class1;eq_class2]==2)=repmat([1,1],s(3),n_eq_idx1);
    scls2(:,n+2,:,[eq_class1;eq_class2]==2)=repmat([1,1],s(3),n_eq_idx2);   
else
    s1=size(rots_grid1);
    s2=size(rots_grid2);
    n_eq_idx1=sum(eq_class1==2);
    n_eq_idx2=sum(eq_class2==2);
    scls1(:,n-1+ceil(n/2),:,eq_class1==2)=repmat([1;1],1,s1(3),n_eq_idx1);
    scls2(:,n-1+ceil(n/2),:,eq_class2==2)=repmat([1;1],1,s2(3),n_eq_idx2);
    n_eq_idx1=sum(eq_class1==3);
    n_eq_idx2=sum(eq_class2==3);
    scls1(:,n-1+ceil(n/2)+1,:,eq_class1==3)=repmat([1;1],1,s1(3),n_eq_idx1);
    scls2(:,n-1+ceil(n/2)+1,:,eq_class2==3)=repmat([1;1],1,s1(3),n_eq_idx2); 
end
%z equators are extremely degenerate. 
% n_eq_idx1=sum(eq_idx==3);
% n_eq_idx2=sum(eq_idx2==3);
% scls1(:,:,:,[eq_class1,eq_class2]==3)=repmat([1,1],2*n-1,s(3),n_eq_idx1);
% scls2(:,:,:,[eq_class1,eq_class2]==3)=repmat([1,1],2*n-1,s(3),n_eq_idx2);

rad=min(round(180/n/4),5); %At least 5 degrees radius
dtheta=round(180*dtheta/pi);
eq_linIdx_lists=cell(ntheta,n_eq1+n_eq2,2); %A cell for each list
for j=nonTv_eq_idx1
    count_eq=count_eq+1;
    for i=1:ntheta
        idx1=circSeq((91+(i-1)*dtheta)-rad,(91+(i-1)*dtheta)+rad,L);
        idx2=circSeq((271+(i-1)*dtheta)-rad,(271+(i-1)*dtheta)+rad,L);
        geq_than_pi_idx=idx2>180;
        idx1(geq_than_pi_idx)=mod(idx1(geq_than_pi_idx)+L/2,L);
        idx2(geq_than_pi_idx)=idx2(geq_than_pi_idx)-L/2;
        idx1(idx1==0)=L;
        idx2(idx2==0)=L;
        eq_linIdx_lists{i,count_eq,1}=sub2ind([L,L/2],idx1,idx2); %DO: uint16 after debug
        eq_linIdx_lists{i,count_eq,2}=idx2; %Since array of eq meaures from allEqMeasures is 1x180
    end
end

count_eq=0;
for j=nonTv_eq_idx2
    count_eq=count_eq+1;
    for i=1:ntheta
        idx1=circSeq((91+(i-1)*dtheta)-rad,(91+(i-1)*dtheta)+rad,L);
        idx2=circSeq((271+(i-1)*dtheta)-rad,(271+(i-1)*dtheta)+rad,L);
        geq_than_pi_idx=idx2>180;
        idx1(geq_than_pi_idx)=mod(idx1(geq_than_pi_idx)+L/2,L);
        idx2(geq_than_pi_idx)=idx2(geq_than_pi_idx)-L/2;
        idx1(idx1==0)=L;
        idx2(idx2==0)=L;
        eq_linIdx_lists{i,count_eq+n_eq1,1}=sub2ind([L,L/2],idx1,idx2); %DO: uint16 after debug
        eq_linIdx_lists{i,count_eq+n_eq1,2}=idx2;
    end
end

%% Compute top view in plane rotation.
% For area top view image, we compute tv measure from the line 180-alpha,
% where alpha is the in plane rotation angle. 
% Though we begin with a theoratical zero in plane rotation for each top
% view image candidate, in practice it's not perfectly top view and so we 
% get a little in plane rotation which we approximate by slcs here. 
% tv_idx1=find(eq_class1>3);
% tv_idx2=find(eq_class2>3);
% n_tv1=length(tv_idx1);
% n_tv2=length(tv_idx2);
% n_tv=sum(eq_class1>3)+sum(eq_class2>3);
% tv_zero_angles=zeros(n_tv,2);
% for k=1:n_tv1
%     switch eq_class1(tv_idx1(k))
%         case 4
%             tv_zero_angles(k,:)=scls1(1,[1,2],1,tv_idx1(k)); %z equator
%         case 5
%             tv_zero_angles(k,:)=scls1(1,[1,3],1,tv_idx1(k)); %y equator
%         case 6
%             tv_zero_angles(k,:)=scls1(1,[2,3],1,tv_idx1(k)); %x equator
%     end
% end
% for k=1:n_tv2
%     switch eq_class2(tv_idx2(k))
%         case 4
%             tv_zero_angles(n_tv1+k,:)=scls2(1,[1,2],1,tv_idx2(k)); %z equator
%         case 5
%             tv_zero_angles(n_tv1+k,:)=scls2(1,[1,3],1,tv_idx2(k)); %y equator
%         case 6
%             tv_zero_angles(n_tv1+k,:)=scls2(1,[2,3],1,tv_idx2(k)); %x equator
%     end
% end
%TO DO: Add this functionality to cryo_clmatrix...

%% Compute linear indices of cls for correlations table
%scls1=scls1*180/pi; %convert to degrees 
s1=size(scls1);
scls1=reshape(scls1,[2,prod(s1(2:end))])';
[scls_lookup1,~]=getClsFromRijs_scls_Dn(scls1);

%scls2=scls2*180/pi; %convert to degrees 
s2=size(scls2);
scls2=reshape(scls2,[2,prod(s2(2:end))])';
[scls_lookup2,~]=getClsFromRijs_scls_Dn(scls2);

%% Top view candidates - A z-axis top view image is an in-plane rotation. 
%we need to register the angle theta of the in-plane rotation. The line
%from which we measure if an image is top view is 180-theta. 
%Similar logic goes for x and y axis top view images. 
% tv_idx=eq_class1>3;
% scls1(1,1,:,tv_idx)=repmat(0:dtheta:(ntheta-1)*dtheta,sum(tv_idx),1)';
% scls1(1,1,:,tv_idx)=mod(pi-scls1(1,1,:,tv_idx),2*pi);
% tv_idx=eq_class2>3;
% scls2(1,1,:,tv_idx)=repmat(0:dtheta:(ntheta-1)*dtheta,sum(tv_idx),1)';
% scls2(1,1,:,tv_idx)=mod(pi-scls2(1,1,:,tv_idx),2*pi);

% scls1=reshape(sub_idx1',s1);
% scls2=reshape(sub_idx2',s2);

% scls1=reshape(scls1,[2,prod(s1(2:end))])';
% scls2=reshape(scls2,[2,prod(s2(2:end))])';
% scls_lookup1=uint16(sub2ind([360,180],scls1(:,1),scls1(:,2)));
% scls_lookup2=uint16(sub2ind([360,180],scls2(:,1),scls2(:,2)));

%% Register some data for self common lines
ntheta=size(rots_grid1,3);
n_eq_idx=sum(eq_idx);
eq_inplane_idx=zeros(ntheta,nrot);
eq_inplane_idx(:,eq_idx)=ones(ntheta,n_eq_idx);
n_eq_idx2=sum(eq_idx2);
eq_inplane_idx2=zeros(ntheta,nrot2);
eq_inplane_idx2(:,eq_idx2)=ones(ntheta,n_eq_idx2);
eq_inplane_idx=eq_inplane_idx(:);
eq_inplane_idx2=eq_inplane_idx2(:);
%tv_idx=[find(eq_class1'>3),nrot+find(eq_class2'>3)];
%n_tv=sum(eq_class1>3)+sum(eq_class2>3);

%Register non equator indices
nNonEq1=length(find(eq_class1'==0));
nNonEq=nNonEq1+length(find(eq_class2'==0));
nonEq_idx=zeros(ntheta,nNonEq);
nonEq_idx(1,:)=([find(eq_class1'==0),nrot+find(eq_class2'==0)]-1)*ntheta+1;
for j=1:ntheta-1
    nonEq_idx(1+j,:)=nonEq_idx(1,:)+j;
end
nonEq_idx=nonEq_idx(:);

scls_lookup_data=struct('scls_lookup1',scls_lookup1,...
    'scls_lookup2',scls_lookup2,'eq_idx',eq_idx,'eq_idx2',eq_idx2,...
    'ntheta',lookup_data.ntheta,'eq_inplane_idx',eq_inplane_idx,...
    'eq_inplane_idx2',eq_inplane_idx2,'eq_class1',eq_class1,...
    'eq_class2',eq_class2,'dtheta',dtheta,'L',L,...
    'nonTv_eq_idx',[nonTv_eq_idx1,nrot+nonTv_eq_idx2],...
    'nonEq_idx',nonEq_idx,'n_eq',n_eq1+n_eq2,...
    'eq2eq_Rij_table11',...
    lookup_data.eq2eq_Rij_table11,'n',n,'eq2eq_Rij_table12',...
    lookup_data.eq2eq_Rij_table12);
scls_lookup_data.eq_idx_lists=eq_linIdx_lists;
%'tv_zero_angles',tv_zero_angles,

