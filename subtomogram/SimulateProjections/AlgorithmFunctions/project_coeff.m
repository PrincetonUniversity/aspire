function [coeff_pos_k] = project_coeff(CR, Crm, InnerP, L, N, ns, nn)

% Project rotated coefficients to simulate images in the real space on a
% uniform Cartesian grid.
%
% INPUT: 
%   CR: rotated spherical expansion coefficients
%   Crm: legendre polynomials Ylm(pi/2,0) or with tilt series adjusted to
%       pi/2
%   InnerP: inner products between spherical Bessel basis in 3D and Fourier
%       Bessel basis in 2D
%   L: size of the original volume and projection
%   N: number of projection images
%   ns: number of slices in tomographic tilt
%   nn: total number of radius basis used across all frequencies.

% OUTPUT:
%   coeff_pos_k: coefficients of the positive Fourier Bessel coefficients
%       after projection

Psilm = cell(1,L); % index l*m, inside cell rotation*q
nls = zeros(L,1);
for ll = 1:L
    nls(ll) = size(InnerP{ll,ll},2);
    Psilm{1, ll} = zeros(N*ns, nls(ll));
end

for ll = 0:L-1
    temp = permute(CR{ll+1,1},[3 1 2]);
    temp = mat2cell(reshape(temp,N,[]),N,size(temp,2)*ones(1,ll+1));
    temp = cellfun(@(x,y,z) bsxfun(@times,x*y,z), temp, InnerP(ll+1,1:ll+1), Crm(ll+1,1:ll+1),'UniformOutput',false);
    temp = cellfun(@(x) reshape(permute(x,[3 1 2]),N*ns,[]), temp,'UniformOutput',false);
    Psilm(1,1:ll+1) = cellfun(@(x,y) x+y, ...
        Psilm(1,1:ll+1), temp, 'UniformOutput', false);
end

coeff_pos_k = zeros(nn,N*ns);
iter = 0;
for ll = 1:L
    nl = nls(ll);
    coeff_pos_k(iter+1:iter+nl,:) = transpose(Psilm{1,ll});
    iter = iter+nl;
end

end

