function [ U_ns, D, ang_freqs, rad_freqs, X ] = FB_SPCA(data, r_max)
%This function computes the steerable PCA basis for the Fourier Bessel
%expansion coefficients
%current code only accept images of odd dimension
%   Input:  data of size LxLxP
%           r_max radius of the mask
%   Output: U_ns: eigenimages
%           D eigenvalues
%           ang_freqs: associated angular frequencies
%           rad_freqs: associated radial frequencies
%           X: mean image
%Zhizhen Zhao Sep 2012

L=size(data, 1);
N=floor(L/2);
P=size(data, 3);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
data=reshape(data, L^2, P);
data=data(r<=r_max, :);

[ Phi_ns, ang_freqs, rad_freqs ]=Bessel_ns_v5(r_max);

%%% Computing Coefficients using least squre
coeff = [Phi_ns, conj(Phi_ns(:, ang_freqs~=0))]\data;
coeff = coeff(1:length(ang_freqs), :);

%%Determine mean of the images:
Mean_Im = real(Phi_ns(:, ang_freqs==0)*mean(coeff(ang_freqs==0, :), 2));
X=zeros(L);
X(r<=r_max)=Mean_Im;

l=size(Phi_ns, 2);
U=zeros(l, l);
D=zeros(l, 1);

for k=1 : max(ang_freqs)+1
    tmp=coeff(ang_freqs==k-1, :);
    if k==1
        mean_tmp=mean(tmp, 2);
        for i=1:P
            tmp(:, i)=tmp(:, i)-mean_tmp;
        end;
    end;
    C=1/(P)*real(tmp*tmp');
    [u, d]=eig(C);
    d=diag(d);
    [d, id]=sort(d, 'descend');
    u=u(:, id);
    id1=find(ang_freqs==k-1, 1, 'first');
    id2=find(ang_freqs==k-1, 1, 'last');
    U(id1:id2, id1:id2)=u;
    D(id1:id2)=d;
end;

D=real(D);
U_ns=Phi_ns*U; 
[D, id]=sort(D, 'descend');
U_ns=U_ns(:, id);

ang_freqs=ang_freqs(id);
rad_freqs=rad_freqs(id);

end
