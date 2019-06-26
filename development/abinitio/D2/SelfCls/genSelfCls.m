
function [scls_lookup_data]=genSelfCls(lookup_data,eqf)


if nargin>1
    eq_filter_angle=eqf;
else
    eq_filter_angle=lookup_data.eq_filter_angle;
end

L=360;
rots_grid1=lookup_data.rots_grid;
rots_grid2=lookup_data.rots_grid2;
sphere_grid1=squeeze(rots_grid1(:,3,:));
sphere_grid2=squeeze(rots_grid2(:,3,:));
[eq_idx,eq_class1]=markEquators(sphere_grid1,eq_filter_angle);
[eq_idx2,eq_class2]=markEquators(sphere_grid2,eq_filter_angle);
rots_grid1=lookup_data.inplane_rotated_grid;
rots_grid2=lookup_data.inplane_rotated_grid2;
dtheta=lookup_data.dtheta;
ntheta=lookup_data.ntheta;

s=size(rots_grid1);
Riis1=zeros([s(1:2),3,s(3:4)]);
scls1=zeros([2,3,s(3:4)]);
scls1_deb=zeros([2,3,360,s(4)]);
nrot=s(4);
for i=1:nrot
    Ris=rots_grid1(:,:,:,i);
    
    gRis=Ris;
    gRis(1:2,:,:)=-gRis(1:2,:,:);
    Riis1(:,:,1,:,i)=multiprod(multi_transpose(Ris),gRis);
    scls1(1,1,:,i)=atan2(Riis1(3,1,1,:,i),-Riis1(3,2,1,:,i));
    scls1(2,1,:,i)=atan2(-Riis1(1,3,1,:,i),Riis1(2,3,1,:,i));
    
    gRis=Ris;
    gRis([1,3],:,:)=-gRis([1,3],:,:);
    Riis1(:,:,2,:,i)=multiprod(multi_transpose(Ris),gRis);
    scls1(1,2,:,i)=atan2(Riis1(3,1,2,:,i),-Riis1(3,2,2,:,i));
    scls1(2,2,:,i)=atan2(-Riis1(1,3,2,:,i),Riis1(2,3,2,:,i));
    
    gRis=Ris;
    gRis(2:3,:,:)=-gRis(2:3,:,:);
    Riis1(:,:,3,:,i)=multiprod(multi_transpose(Ris),gRis);
    scls1(1,3,:,i)=atan2(Riis1(3,1,3,:,i),-Riis1(3,2,3,:,i));
    scls1(2,3,:,i)=atan2(-Riis1(1,3,3,:,i),Riis1(2,3,3,:,i));
    
    scls1_deb(1,1,:,i)=scls1(1,1,1,i)+(0:359);
    scls1_deb(2,1,:,i)=scls1(2,1,1,i)+(0:359);
    scls1_deb(1,2,:,i)=scls1(1,2,1,i)+(0:359);
    scls1_deb(2,2,:,i)=scls1(2,2,1,i)+(0:359);
    scls1_deb(1,3,:,i)=scls1(1,3,1,i)+(0:359);
    scls1_deb(2,3,:,i)=scls1(2,3,1,i)+(0:359);
    %scls1_deb(:,:,:,i)=mod(scls1_deb(:,:,:,i),2*pi);
end

s=size(rots_grid2);
Riis2=zeros([s(1:2),3,s(3:4)]);
scls2=zeros([2,3,s(3:4)]);
scls2_deb=zeros([2,3,360,s(4)]);
nrot2=s(4);
for i=1:nrot2
    Ris=rots_grid2(:,:,:,i);
    
    gRis=Ris;
    gRis(1:2,:,:)=-gRis(1:2,:,:);
    Riis2(:,:,1,:,i)=multiprod(multi_transpose(Ris),gRis);
    scls2(1,1,:,i)=atan2(Riis2(3,1,1,:,i),-Riis2(3,2,1,:,i));
    scls2(2,1,:,i)=atan2(-Riis2(1,3,1,:,i),Riis2(2,3,1,:,i));
    
    gRis=Ris;
    gRis([1,3],:,:)=-gRis([1,3],:,:);
    Riis2(:,:,2,:,i)=multiprod(multi_transpose(Ris),gRis);
    scls2(1,2,:,i)=atan2(Riis2(3,1,2,:,i),-Riis2(3,2,2,:,i));
    scls2(2,2,:,i)=atan2(-Riis2(1,3,2,:,i),Riis2(2,3,2,:,i));
    
    gRis=Ris;
    gRis(2:3,:,:)=-gRis(2:3,:,:);
    Riis2(:,:,3,:,i)=multiprod(multi_transpose(Ris),gRis);
    scls2(1,3,:,i)=atan2(Riis2(3,1,3,:,i),-Riis2(3,2,3,:,i));
    scls2(2,3,:,i)=atan2(-Riis2(1,3,3,:,i),Riis2(2,3,3,:,i));
    
    scls2_deb(1,1,:,i)=scls2(1,1,1,i)+(0:359);
    scls2_deb(2,1,:,i)=scls2(2,1,1,i)+(0:359);
    scls2_deb(1,2,:,i)=scls2(1,2,1,i)+(0:359);
    scls2_deb(2,2,:,i)=scls2(2,2,1,i)+(0:359);
    scls2_deb(1,3,:,i)=scls2(1,3,1,i)+(0:359);
    scls2_deb(2,3,:,i)=scls2(2,3,1,i)+(0:359);
end

%% Prepare common line coordinates
scls1=mod(scls1,2*pi); %make all angles non negative 
scls2=mod(scls2,2*pi); %make all angles non negative

%  Deal with non top view equators
%  put 'real' self common line at 2 first coordinates, candidate for
%  perpendicular line is in 3rd coordinate. 
%  The real 2 common lines of an equator are very close but anti-podal so
%  we flip one of them to make them nearly the same
%  We try all possible line candidates between the 2 'real' common lines

nonTv_xEq_idx=(eq_class1==1);
scls1(:,:,:,nonTv_xEq_idx)=scls1(:,[2,3,1],:,nonTv_xEq_idx);
scls1(:,1,:,nonTv_xEq_idx)=scls1([2,1],1,:,nonTv_xEq_idx);

nonTv_yEq_idx=(eq_class1==2);
scls1(:,:,:,nonTv_yEq_idx)=scls1(:,[1,3,2],:,nonTv_yEq_idx);
scls1(:,1,:,nonTv_yEq_idx)=scls1([2,1],1,:,nonTv_yEq_idx);

%No need to rearrange z but flip common line to antipodals
nonTv_zEq_idx=(eq_class1==3);
scls1(:,1,:,nonTv_zEq_idx)=scls1([2,1],1,:,nonTv_zEq_idx);

nonTv_xEq_idx=(eq_class2==1);
scls2(:,:,:,nonTv_xEq_idx)=scls2(:,[2,3,1],:,nonTv_xEq_idx);
scls2(:,1,:,nonTv_xEq_idx)=scls2([2,1],1,:,nonTv_xEq_idx);

nonTv_yEq_idx=(eq_class2==2);
scls2(:,:,:,nonTv_yEq_idx)=scls2(:,[1,3,2],:,nonTv_yEq_idx);
scls2(:,1,:,nonTv_yEq_idx)=scls2([2,1],1,:,nonTv_yEq_idx);

%No need to rearrange z but flip common line to antipodals
nonTv_zEq_idx=(eq_class2==3);
scls2(:,1,:,nonTv_zEq_idx)=scls2([2,1],1,:,nonTv_zEq_idx);

%Make sure angle range is <= 180 degrees, first scls1
nonTv_eq_idx=eq_class1>0 & eq_class1<4;
p1=scls1(:,1,:,nonTv_eq_idx)>scls1(:,2,:,nonTv_eq_idx);
p1=p1(1,1,:,:) & p1(2,1,:,:);
p2=scls1(:,1,:,nonTv_eq_idx)-scls1(:,2,:,nonTv_eq_idx)<-pi;
p2=p2(1,1,:,:) | p2(2,1,:,:);
p=p1 | p2;
flipped_scls1=scls1(:,[1,2],:,nonTv_eq_idx);
p=repmat(p,2,2,1,1);
flipped_scls1(:,[1,2],:,:)=p.*scls1(:,[2,1],:,nonTv_eq_idx);
scls1(:,[1,2],:,nonTv_eq_idx)=~p.*scls1(:,[1,2],:,nonTv_eq_idx)+flipped_scls1;

%Now for scls2 
nonTv_eq_idx=eq_class2>0 & eq_class2<4;
p1=scls2(:,1,:,nonTv_eq_idx)>scls2(:,2,:,nonTv_eq_idx);
p1=p1(1,1,:,:) & p1(2,1,:,:);
p2=scls2(:,1,:,nonTv_eq_idx)-scls2(:,2,:,nonTv_eq_idx)<-pi;
p2=p2(1,1,:,:) | p2(2,1,:,:);
p=p1 | p2;
flipped_scls2=scls2(:,[1,2],:,nonTv_eq_idx);
p=repmat(p,2,2,1,1);
flipped_scls2(:,[1,2],:,:)=p.*scls2(:,[2,1],:,nonTv_eq_idx);
scls2(:,[1,2],:,nonTv_eq_idx)=~p.*scls2(:,[1,2],:,nonTv_eq_idx)+flipped_scls2;

%Finally transform angles from radians to degrees. 
scls1=round(scls1*180/pi)+1;
scls2=round(scls2*180/pi)+1;

%% Create candidate common line linear indices lists for equators
nonTv_eq_idx1=find(eq_class1>0 & eq_class1<4)';
n_eq1=length(nonTv_eq_idx1);
nonTv_eq_idx2=find(eq_class2>0 & eq_class2<4)';
n_eq2=length(nonTv_eq_idx2);
count_eq=0;

eq_linIdx_lists=cell(ntheta,n_eq1+n_eq2,2); %A cell for each list
for j=nonTv_eq_idx1
    count_eq=count_eq+1;
    for i=1:ntheta
        idx1=circSeq(scls1(1,1,i,j),scls1(1,2,i,j),L);
        idx2=circSeq(scls1(2,1,i,j),scls1(2,2,i,j),L);
        geq_than_pi_idx=idx2>180;
        idx1(geq_than_pi_idx)=mod(idx1(geq_than_pi_idx)+L/2,L);
        idx2(geq_than_pi_idx)=idx2(geq_than_pi_idx)-L/2;
        idx1(idx1==0)=L;
        idx2(idx2==0)=L;
        eq_linIdx_lists{i,count_eq,1}=sub2ind([L,L/2],idx1,idx2); %DO: uint16 after debug
        eq_linIdx_lists{i,count_eq,2}=idx2; %Since array of eq meaures from allEqMeasures is 1x180
    end
end

%DEBUG
% scls1=scls1_deb;
% scls2=scls2_deb;
%END DEBUG
count_eq=0;

for j=nonTv_eq_idx2
    count_eq=count_eq+1;
    for i=1:ntheta
        idx1=circSeq(scls2(1,1,i,j),scls2(1,2,i,j),L);
        idx2=circSeq(scls2(2,1,i,j),scls2(2,2,i,j),L);
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
tv_idx1=find(eq_class1>3);
tv_idx2=find(eq_class2>3);
n_tv1=length(tv_idx1);
n_tv2=length(tv_idx2);
n_tv=sum(eq_class1>3)+sum(eq_class2>3);
tv_zero_angles=zeros(n_tv,2);
for k=1:n_tv1
    switch eq_class1(tv_idx1(k))
        case 4
            tv_zero_angles(k,:)=scls1(1,[1,2],1,tv_idx1(k)); %z equator
        case 5
            tv_zero_angles(k,:)=scls1(1,[1,3],1,tv_idx1(k)); %y equator
        case 6
            tv_zero_angles(k,:)=scls1(1,[2,3],1,tv_idx1(k)); %x equator
    end
end
for k=1:n_tv2
    switch eq_class2(tv_idx2(k))
        case 4
            tv_zero_angles(n_tv1+k,:)=scls2(1,[1,2],1,tv_idx2(k)); %z equator
        case 5
            tv_zero_angles(n_tv1+k,:)=scls2(1,[1,3],1,tv_idx2(k)); %y equator
        case 6
            tv_zero_angles(n_tv1+k,:)=scls2(1,[2,3],1,tv_idx2(k)); %x equator
    end
end
%TO DO: Add this functionality to cryo_clmatrix...

%% Compute linear indices of cls for correlations table
%scls1=scls1*180/pi; %convert to degrees 
s1=size(scls1);
scls1=reshape(scls1,[2,prod(s1(2:end))])';
[scls_lookup1,~]=getClsFromRijs2_scls(scls1);

%scls2=scls2*180/pi; %convert to degrees 
s2=size(scls2);
scls2=reshape(scls2,[2,prod(s2(2:end))])';
[scls_lookup2,~]=getClsFromRijs2_scls(scls2);

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
tv_idx=[find(eq_class1'>3),nrot+find(eq_class2'>3)];
n_tv=sum(eq_class1>3)+sum(eq_class2>3);

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
    'nonEq_idx',nonEq_idx,'tv_idx',tv_idx,'n_eq',n_eq1+n_eq2,'n_tv',n_tv,...
    'tv_zero_angles',tv_zero_angles,'eq2eq_Rij_table11',...
    lookup_data.eq2eq_Rij_table11,'eq2eq_Rij_table12',...
    lookup_data.eq2eq_Rij_table12);
scls_lookup_data.eq_idx_lists=eq_linIdx_lists;

