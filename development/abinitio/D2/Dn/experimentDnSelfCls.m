
n=5;
theta=3*pi/n;
v=[-sin(theta);cos(theta);0]; v=v/norm(v);
u=[0;0;1];
a=1;b=1.5;
w=b*v+a*u; w=w/norm(w);
projDir=[cos(theta);sin(theta);0];
R=[cross(projDir,w)/norm(cross(projDir,w)),projDir,w];

%Z equator
theta=10*pi/180;
projDir=[cos(theta);sin(theta);0];
v=[-sin(theta);cos(theta);0]; v=v/norm(v);
S=[cross(v,projDir)/norm(cross(v,projDir)),v,projDir];

[gs,g_z]=gns(n);
gs=cat(3,g_z,gs);
dirs=zeros(3,n+1);
for k=1:n+1
    dirs(:,k)=gs(:,:,k)*w;
end


RQ=zeros(3,3,n+1);
cls=zeros(n+1,2);
for k=1:n+1
    RQ(:,:,k)=R'*gs(:,:,k)*R;
    f=1;%sqrt(1-RQ(3,3)^2);
    cls(k,2)=atan2(-RQ(3,1,k)/f,RQ(3,2,k)/f);
    cls(k,1)=atan2(RQ(1,3,k)/f,-RQ(2,3,k)/f);
end

projR = cryo_project(vol,R');
projR = projR';
[pf,~]=cryo_pft(projR,size(projR,1),360);
pf=cat(3,pf,pf);
[C,est_shifts]=cryo_calc_proj_corrs(pf,0,0,1);
[cls2,corrs]=cryo_plot_cls_Dn(R,R,5,C);

[g1,g2]=gns(n);
RQ=zeros(3,3,2*n);
cls=zeros(2*n,2);
for k=1:n
    RQ(:,:,k)=R'*g1(:,:,k)*R;
    RQ(:,:,k+n)=R'*g2(:,:,k)*R;
    f=1;%sqrt(1-RQ(3,3)^2);
    cls(k,2)=atan2(-RQ(3,1,k)/f,RQ(3,2,k)/f);
    cls(k,1)=atan2(RQ(1,3,k)/f,-RQ(2,3,k)/f);
    cls(k+n,2)=atan2(-RQ(3,1,k+n)/f,RQ(3,2,k+n)/f);
    cls(k+n,1)=atan2(RQ(1,3,k+n)/f,-RQ(2,3,k+n)/f);
end
cls=round(cls*180/pi);
cls=mod(cls+180,360);

%Compute self common line impersonator for D_{2n-1}
eqProjDir=R(:,3);
u=cross(eqProjDir,[0;0;1]);
u=u/norm(u);
if u(2)<0
    u=-u;
end
scl_p=cross(eqProjDir,u);
scl_p=scl_p/norm(scl_p);
tmp_c=R'*scl_p;
scl_p_c=atan2(tmp_c(2),tmp_c(1));

%Compute eq_measure
tmp2=tmp(:,11:100)-flip(tmp(:,102:191),2);

grid_in=cat(3,R,Q);
[projs,Rijs_gt,q,ref_shifts]=genDataForSimulation(vol,...
    100,max_shift,1,1000000,s,0,[],grid_in);


