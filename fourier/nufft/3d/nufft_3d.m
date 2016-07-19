function f=nufft_3d(alpha,omega,precision,M)
%
% Optimized version of nufft_3d_ref.m.
% See nufft_3d_ref.m. for more information.
%
% Yoel Shkolnisky, January 2010.

if nargin<3
    precision='single';
end

[b,m,q]=nufftparams(precision);

if length(size(omega))~=2
    error('omega must be a 3 column array')
end

if size(omega,2)~=3
    error('omega must be a 3 column array');
end

alpha=alpha(:);
n=numel(alpha);
mu=round(omega.*m);

if size(omega,1)~=n
    error('alpha and omega must have the same number of rows');
end

% All precomputations are borken into 1D ones.
Px=zeros(n,q+1);
Py=zeros(n,q+1);
Pz=zeros(n,q+1);
for j=-q/2:q/2
    tmp1=(m.*omega(:,1)-(mu(:,1)+j)).^2;
    tmp2=((m.*omega(:,2)-(mu(:,2)+j)).^2);
    tmp3=((m.*omega(:,3)-(mu(:,3)+j)).^2);
    Px(:,j+q/2+1)=exp(-tmp1/(4*b))./(2*sqrt(b*pi));
    Py(:,j+q/2+1)=exp(-tmp2/(4*b))./(2*sqrt(b*pi));
    Pz(:,j+q/2+1)=exp(-tmp3/(4*b))./(2*sqrt(b*pi));
end

tau=nufft3dauxmx(n,M,m,q,mu,Px.',Py.',Pz.',alpha);

tau = ifftshift(ifftshift(ifftshift(tau, 1), 2), 3);
T = ifftn(tau);
T = fftshift(fftshift(fftshift(T, 1), 2), 3);
T=T*numel(T);

low_idx_M=-ceil((M-1)/2);
high_idx_M=floor((M-1)/2);
idx=low_idx_M:high_idx_M;
E=exp(b.*(2.*pi.*idx./(m*M)).^2);
E=E(:);

offset1=ceil((m*M-1)/2);
offset2=offset1+low_idx_M+1;
E2=E*E.';
E3=reshape(E2(:)*E.',[M M M]);
f=T(offset2:offset2+M-1,offset2:offset2+M-1,offset2:offset2+M-1).*E3;
