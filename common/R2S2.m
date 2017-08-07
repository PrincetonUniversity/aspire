function PR=R2S2(Rs,Ntheta)
%
% Compute direction vectors of Fourier rays correspoding to given
% rotation matrices.
%
% R is a 3x3xK array containing K rotation matrices.
%
% For each rotation matrix, which corresponds to the projection
% orientation of one of the projections, the function computes the
% direction vectors of the Fourier rays in that projection. 
% There are Ntheta (default 300) Fourier rays in each projection. The size
% of the output is thus K*Nthetax3.
% 
% Examples:
%   rotations = rand_rots(100);
%   L = 36;
%   dir = R2S2(rotations, L);
%
% Yoel Shkolnisky, August 2010.

if nargin<2
    Ntheta=300;
end

n=size(Rs,3);
dtheta=2*pi/Ntheta;

PR=zeros(Ntheta*n,3);
for k=1:n
    R=Rs(:,:,k);
    e1=R*([1 0 0].');   % We don't really need to multiply - just for
    e2=R*([0 1 0].');   % explanatory purposes.

    for j=1:Ntheta
        PR((k-1)*Ntheta+j,:)=cos((j-1)*dtheta).*e1+sin((j-1)*dtheta).*e2;
    end

end

return;
