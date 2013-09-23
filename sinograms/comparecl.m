function [ prob, clean_clmatrix, n_cl ] = comparecl( clstack, ...
    ref_clmatrix, n_theta, max_angle )
% compare the estimated clstack to the true ref_clmatrix and compute the
% correctness.
%
% Input: 
%   clstack: estimated commonline stack.
%   ref_clmatrix: true commonline stack.
%   n_theta: angular resolution.
%   max_angle: maximum tolerant error
%
% Output: 
%   prob: the probability of good commonlines
%   clean_clmatrix: the commonline stack with good commonlines
%   n_cl: the number of good commonlines for each image.
%
% Lanhui Wang, July 3, 2013 (Based on cryo_clmatrix_v3.m by Yoel)


PI=4*atan(1.0);
max_angle=max_angle*PI/180;
n_theta=n_theta/2;
angle_tol=2*sin(max_angle/2)+1.0e-10; % Add 1.0e-10 to account for roundoff error
alpha=2*PI*sqrt(-1)/(2*n_theta);
n_proj=size(clstack,1);
correct=0;
clean_clmatrix=zeros(size(clstack));
sum_cl=0;

for k1=1:n_proj
    for k2=k1+1:n_proj
    tcl1=ref_clmatrix(k1,k2);
    tcl2=ref_clmatrix(k2,k1);
    if tcl1~=0
        sum_cl=sum_cl+1;
    l1=clstack(k1,k2);
    l2=clstack(k2,k1);
    d1s=abs(exp(alpha*(l1-1))-exp(alpha*(tcl1-1)));
    d2s=abs(exp(alpha*(l2-1))-exp(alpha*(tcl2-1)));
    d1f=abs(exp(alpha*(l1-1)+sqrt(-1)*PI)-exp(alpha*(tcl1-1)));
    d2f=abs(exp(alpha*(l2-1)+sqrt(-1)*PI)-exp(alpha*(tcl2-1)));

    if (d1s<=angle_tol) && (d2s<=angle_tol) || ...
        (d1f<=angle_tol) && (d2f<=angle_tol)
        correct=correct+1;
        clean_clmatrix(k1,k2)=l1;
        clean_clmatrix(k2,k1)=l2;
    end
    end
    end
end
prob=correct/sum_cl;

n_cl=zeros(1,n_proj);
for k=1:n_proj
    n_cl(k)=length(find(clean_clmatrix(k,:)));
end
end

