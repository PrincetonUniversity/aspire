function f=nufft_3d_ref(alpha,omega,precision,M)
%
% Non-equally spaced FFT.
%
% The function computes the sums
%                n
%    f(j) = sum alpha(k)*exp(2*pi*i*(dot(j,omega(k))/M)
%               k=1
% for j in [-M/2,...,M/2]^3, where
% omega(k)=(omega1(k),omega2(k),omega3(k)), with omega(k,l) in [-n/2,n/2].
%
% The complexity of the algorithm is O(n*log(1/eps)+(m*M)^3(log(M))).
% m is the oversampling factor (m=2).
%
% Input parameters:
%    alpha    Coefficients in the sums above. Real or complex numbers.
%             A vector of length n. 
%    omega    Sampling frequnecies. Real numbers in the range [-n/2,n/2]. 
%             Must be an array with three columns.
%    precision  'double' or 'single'. Default is 'single'.
%
% omega(k) is the frequency that corresponds to the Fourier coefficient
% alpha(k). alpha and omega must have the same length. If alpha is given as
% a three-dimensional array with n elements, then omega(k,1) is in the
% direction of the fast dimension (rows for a matrix) , omega(k,2) is
% second fastest dimension, and so on.
%
% The code was tested for M<n and everything seems to work fine.
%
% Note: 
% This code is a rather slow implementation of the NUFFT algorithm. It
% follows exactly the description in the paper. For better performance,
% call nufft_3d.m.
%
% Output:
%    f        The sums defined above.
%
% Yoel Shkolnisky, December 2009.
%

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

tau=zeros(m*M,m*M,m*M);
offset1=ceil((m*M-1)/2);
for j1=-q/2:q/2
    for j2=-q/2:q/2
        for j3=-q/2:q/2
            idx=mu+repmat([j1 j2 j3],n,1);
            idx=mod(idx+offset1,m*M)+1;
            idx=(idx(:,3)-1)*(m*M)^2+(idx(:,2)-1)*m*M+idx(:,1); % convert to linear index
            tau(idx(:))=tau(idx(:))+...
                Px(:,j1+q/2+1).*Py(:,j2+q/2+1).*Pz(:,j3+q/2+1).*alpha(:);
        end
    end
end

T=fftshift(ifftn(ifftshift(tau)));
T=T*numel(T);

low_idx_M=-ceil((M-1)/2);
high_idx_M=floor((M-1)/2);
idx=low_idx_M:high_idx_M;
E=exp(b.*(2.*pi.*idx./(m*M)).^2);
E=E(:);

offset2=offset1+low_idx_M+1;

f=zeros(M,M,M);
for j1=1:M
    for j2=1:M
        for j3=1:M
            f(j1,j2,j3)=T(offset2+j1-1,offset2+j2-1,offset2+j3-1).*E(j1).*E(j2)*E(j3);
        end
    end
end
