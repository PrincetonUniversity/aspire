
function [colUnmix,best_score,best_t]=unmixColors_best(colors,w_in)
%Unmix vectors
v=colors(:,1);
l=length(v);
v=reshape(v,3,l/3);
norms=sqrt(sum(v.^2,1));
v=v./norms;

x=colors(:,2);
x=reshape(x,3,l/3);
norms=sqrt(sum(x.^2,1));
x=x./norms;

p=perms(1:3);
u=[1;0;-1];
u=u/norm(u);
w=w_in;%[1;-2;1];
w=w/norm(w);
up=zeros(3,6);
wp=up;
for i=1:6
    up(:,i)=u(p(i,:));
    wp(:,i)=w(p(i,:));
    %up(:,i+6)=-u(p(i,:));
    %wp(:,i+6)=-w(p(i,:));
end
up=repmat(up,l/3,1);
wp=repmat(wp,l/3,1);
up=reshape(up,3,l/3,6);
wp=reshape(wp,3,l/3,6);
%up=reshape(up,3,l/3,12);
%wp=reshape(wp,3,l/3,12);

K=[0 -1/sqrt(3) 1/sqrt(3);
   1/sqrt(3) 0 -1/sqrt(3);
   -1/sqrt(3) 1/sqrt(3) 0];
R_theta=@(t)(eye(3)+sin(t)*K+(1-cos(t))*K*K);
best_t=0;
best_score=inf;
idx=0;
ones_l=ones(l,1);
for t=0.5:0.5:360
    idx=idx+1;
    norms=zeros(24,l/3);
    R_t1=R_theta(pi*t/180);
    R_t2=R_theta(pi*(360-t)/180);
    v_rot1=R_t1*v;
    v_rot2=R_t2*v;
    
    x_rot1=R_t1*x;
    x_rot2=R_t2*x;    
    
    for k=1:6
        norms(k,:)=sum((v_rot1-up(:,:,k)).^2,1)+...
            sum((x_rot1+wp(:,:,k)).^2,1);
        norms(k+6,:)=sum((v_rot1-up(:,:,k)).^2,1)+...
            sum((x_rot2+wp(:,:,k)).^2,1);
        
        norms(k+12,:)=sum((v_rot2-up(:,:,k)).^2,1)+...
            sum((x_rot2+wp(:,:,k)).^2,1);  
        norms(k+18,:)=sum((v_rot2-up(:,:,k)).^2,1)+...
            sum((x_rot1+wp(:,:,k)).^2,1);
        
        
        
%         norms(k+12,:)=sum((v_rot1-up(:,:,k)).^2,1)+...
%             sum((x_rot1+wp(:,:,k)).^2,1);
%         norms(k+18,:)=sum((v_rot2-up(:,:,k)).^2,1)+...
%             sum((x_rot2+wp(:,:,k)).^2,1);
    end
    [norms,I]=min(norms,[],1);        
    
    score_t=sum(norms);
    if score_t<best_score
        best_score=score_t;
        best_max=max(norms);
        norms_out=norms;
        best_t=t;
        %colUnmix=v_rot1(:);
        whichTheta=I<13;
        whichTheta=repmat(whichTheta,3,1);
        whichTheta=whichTheta(:);
        colUnmix=(whichTheta.*v_rot1(:)+(ones_l-whichTheta).*v_rot2(:));
    end
end
colUnmix=colUnmix/norm(colUnmix);



