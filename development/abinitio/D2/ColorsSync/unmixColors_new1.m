
function [colors,best_unmix]=unmixColors_new1(ev,ntheta)

np=size(ev,1)/3;
dtheta=360/ntheta;
maxt=360/dtheta+1;
R_theta=@(theta)([cos(theta),-sin(theta);sin(theta),cos(theta)]);
s=inf;
scores=zeros(1,maxt);
idx=0;
for t=0:0.5:maxt
    idx=idx+1;
    unmix_ev=ev*R_theta(pi*t/180);
    s1=reshape(unmix_ev(:,1),3,np);
    [s1,p1]=sort(s1,1,'descend');
    s2=abs(reshape(unmix_ev(:,2),3,np));
    [s2,p2]=sort(s2,1,'descend');
    score11=sum((s1(1,:)+s1(3,:)).^2+s1(2,:).^2);
    score12=sum((2*s2(1,:)-s2(3,:)).^2+(2*s2(2,:)-s2(3,:)).^2+...
        (s2(1,:)-s2(3,:)).^2);
    score21=sum((s2(1,:)+s2(3,:)).^2+s2(2,:).^2);
    score22=sum((2*s1(1,:)-s1(3,:)).^2+(2*s1(2,:)-s1(3,:)).^2+...
        (s1(1,:)-s1(3,:)).^2);
    [scores(idx),whichVec]=min([score11+score12,score21+score22]);
    if scores(idx)<s
        s=scores(idx);
        if whichVec==1
            p=p1;
        else
            p=p2;
        end
        best_unmix=unmix_ev(:,whichVec);
        %maxt=t;
    end    
end

%Assign integers between 1:3 to permutations
colors=zeros(3,np);
for i=1:np
    p_i=p(:,i);
    p_i_sqr=p_i(p_i);
    if sum((p_i_sqr-[1;2;3]).^2)==0 % non cyclic pemutation
        colors(:,i)=p_i;
    else
        colors(:,i)=p_i_sqr;
    end
end
colors=colors(:);









