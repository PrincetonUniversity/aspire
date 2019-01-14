
%Input:  c1,c2=angles in degrees of commons line on images 1 and 2
%        L = Sampling resolution        
%        C = LxL/2 fourier ray correlations matrix
function m=measureEq2(C,c1,L)

c3=mod(c1+L/4,L);%round((c1+c2)/2);
if c3<91
    c3=c3+180;
end
idx=mod([(c3-(1:90))',(c3+(1:90))'+L/2],360);
zero_idx1=abs(idx(:,1))<1e-7;
zero_idx2=abs(idx(:,2))<1e-7;
idx(zero_idx1,1)=360;
idx(zero_idx2,2)=360;

bigger_than_180=idx(:,2)>180;
idx(bigger_than_180,1)=mod(idx(bigger_than_180,1)+180,360);
zero_idx1=abs(idx(:,1))<1e-7;
idx(zero_idx1,1)=360;
idx(bigger_than_180,2)=idx(bigger_than_180,2)-180;
corrs_idx=sub2ind([L,L/2],idx(:,1),idx(:,2));

m=zeros(1,3);
m(1)=gather(median(C(corrs_idx))); %Consider also mean
m(2)=gather(min(C(corrs_idx)));

if c1<91
    c3=c1+180;
end
normal2scl_idx=mod([(c3-(85:95))',(c3-(85:95))'+L/2],360);
zero_idx1=abs(normal2scl_idx(:,1))<1e-7;
zero_idx2=abs(normal2scl_idx(:,2))<1e-7;
normal2scl_idx(zero_idx1,1)=360;
normal2scl_idx(zero_idx2,2)=360;

bigger_than_180=normal2scl_idx(:,2)>180;
normal2scl_idx(bigger_than_180,1)=mod(normal2scl_idx(bigger_than_180,1)+180,360);
zero_idx1=abs(normal2scl_idx(:,1))<1e-7;
normal2scl_idx(zero_idx1,1)=360;
normal2scl_idx(bigger_than_180,2)=normal2scl_idx(bigger_than_180,2)-180;
corrs_idx=sub2ind([L,L/2],normal2scl_idx(:,1),normal2scl_idx(:,2));
m(3)=gather(max(C(corrs_idx)));
    


