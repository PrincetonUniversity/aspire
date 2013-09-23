function f=nufft_2d(alpha,omega,precision,M)
%
% Non-equally spaced FFT.
%
% The function computes the sums
%                n
%    f(j1,j2) = sum alpha(k)*exp(2*pi*i*(j1*omega1(k)+j2*omega2(k))/M)
%               k=1
% for j1,j2=-M/2,...,M/2, where omega(k)=(omega1(k),omega2(k)), with
% omega1(k) and omega2(k) in [-n/2,n/2].
%
% The complexity of the algorithm is O(n*log(1/eps)+m*M(log(M))).
% m is the oversampling factor (m=2).
%
% Input parameters:
%    alpha    Coefficients in the sums above. Real or complex numbers.
%             A vector of length n. 
%    omega    Sampling frequnecies. Real numbers in the range [-n/2,n/2]. 
%             Must be an array with two columns.
%    precision  'double' or 'single'. Default is 'single'.
%
% omega(k) is the frequency that corresponds to the Fourier coefficient
% alpha(k). alpha and omega must have the same length. If alpha is given as
% a two-dimensional array with n elements, then omega(k,1) is in the
% direction of the rows and omega(k,2) is in the direction of the colummns.
%
% The code was tested for M<n and everything seems to work fine.
% The speed of the function can be significantly improved by MEXing the
% gridding loop. See the 3D version nufft3dauxmx.cpp.
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
    error('omega must be a 2 column array')
end

if size(omega,2)~=2
    error('omega must be a 2 column array');
end

alpha=alpha(:);
n=numel(alpha);

if size(omega,1)~=n
    error('alpha and omega must have the same number of rows');
end

mu=round(omega.*m);

P=zeros(n,q+1,q+1);
for j1=-q/2:q/2
    for j2=-q/2:q/2
        tmp=-((m.*omega(:,1)-(mu(:,1)+j1)).^2+(m.*omega(:,2)-(mu(:,2)+j2)).^2)/(4*b);
        P(:,j1+q/2+1,j2+q/2+1)=exp(tmp)./(4*b*pi);
    end
end

% There is no reason to keep P as a 3D array. Everything can be implemeted
% using 1D interpolations.
Px=zeros(n,q+1);
Py=zeros(n,q+1);
for j=-q/2:q/2
    tmp1=(m.*omega(:,1)-(mu(:,1)+j)).^2;
    tmp2=(m.*omega(:,2)-(mu(:,2)+j)).^2;
    Px(:,j+q/2+1)=exp(-tmp1/(4*b))./(2*sqrt(b*pi));
    Py(:,j+q/2+1)=exp(-tmp2/(4*b))./(2*sqrt(b*pi));
end

tau=zeros(m*M,m*M,1);
offset1=ceil((m*M-1)/2);
% for k=1:n
%     for j1=-q/2:q/2
%         for j2=-q/2:q/2
%             idx=mu(k,:)+[j1 j2];
%             idx=mod(idx+offset1,m*M)+1;
%             tau(idx(1),idx(2))=tau(idx(1),idx(2))+P(k,j1+q/2+1,j2+q/2+1)*alpha(k);
%         end
%     end
% end

% Gridding loop start
for j1=-q/2:q/2
    idx1=mu(:,1)+j1;
    idx1=mod(idx1+offset1,m*M)+1;
    W=Px(:,j1+q/2+1).*alpha(:);
    for j2=-q/2:q/2
        idx2=mu(:,2)+j2;
        idx2=mod(idx2+offset1,m*M)+1;
        idx=(idx2-1)*m*M+idx1; % Convert to linear index.
        for k=1:n
            tau(idx(k))=tau(idx(k))+Py(k,j2+q/2+1).*W(k);
        end
    end
end
% Gridding loop end

T=fftshift(ifft2(ifftshift(tau)));
T=T*numel(T);

low_idx_M=-ceil((M-1)/2);
high_idx_M=floor((M-1)/2);
idx=low_idx_M:high_idx_M;
E=exp(b.*(2.*pi.*idx./(m*M)).^2);
E=E(:);
E=E*E.';

offset2=offset1+low_idx_M+1;
f=T(offset2:offset2+M-1,offset2:offset2+M-1).*E;

return;
