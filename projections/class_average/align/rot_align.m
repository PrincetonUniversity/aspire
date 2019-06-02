function [ corr, rot, ITER ]=rot_align(m, coeff, list)
%   rotationally align pairs of images
%   Input:
%       m: frequencies computed from steerable basis
%       coeff: in cell structure, each cell contains coefficients with one
%       angular frequencies. 
%       list: of size Px2 where P is the number of pairs of images. 
%   Output:
%       rot: of size Px1 is the rotatinal alignment angle in degrees
%       ITER: the number of iterations
%   Zhizhen Zhao Feb 10 2012
N_theta=360;
P=size(list, 1);
max_m=max(m);
C=zeros(max_m, P);
m_list=[1:max_m]';
ITER=zeros(P, 1);

for i=1:max_m+1
    C(i, :)=sum(conj(coeff{i}(:, list(1:P, 1))).*coeff{i}(:, list(1:P, 2)), 1);
end
C2=flipud(conj(C(2:end, :)));
B=real((2*max_m+1)*icfftd([C2; C], 1));
[~, rot]=max(B, [], 1);
rot=(rot-max_m-1)*N_theta/(2*max_m+1);

%%Newton-Raphson method for root finding
f_prime=@(x) sum(real(sqrt(-1)*repmat(m_list*(pi/180), 1, P).*C(2:end, P).*exp(sqrt(-1)*m_list*x*pi/180)));
f_prime2=@(x) sum(real((-1)*repmat(m_list.^2*(pi/180)^2, 1, P).*C(2:end, P).*exp(sqrt(-1)*m_list*x*pi/180)));
x_old=-ones(1, P);
x_new=rot;
precision=1e-12;
iter=0;
while max(abs(x_new - x_old)) > precision
    max(abs(x_new - x_old))
    x_old = x_new;
    delta=f_prime(x_old)./f_prime2(x_old);
    x_new = x_new - delta;
    iter=iter+1;
    if iter>10, break, end
end
rot=x_new;
m_list=[0; m_list];
C=C.*exp(sqrt(-1)*m_list*rot*pi/180);
corr=(real(C(1, :))+2*sum(real(C(2:end, :))))/2;

end
