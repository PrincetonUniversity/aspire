function PR=Q2S2(Q,Ntheta)
%
% Compute direction vectors of Fourier rays correspoding to given
% quaternions.
%
% Q is a 4xK array containing quaternions.
%
% For each quaternion, which corresponds to the projection
% orientation of one of the projections, the function computes the
% direction vectors of the Fourier rays in that projection. 
% There are Ntheta (default 300) Fourier rays in each projection. The size
% of the output is thus N*Nthetax3.
% 
% Examples:
%   directions=Q2S2(Q);
%
% Yoel Shkolnisky, September 2008.

if nargin<2
    Ntheta=300;
end

n=size(Q,2);
dtheta=2*pi/Ntheta;

PR=zeros(Ntheta*n,3);
for k=1:n
    R=q_to_rot(Q(:,k));
    e1=inv(R)*([1 0 0].');   % We don't really need to multiply - just for
    e2=inv(R)*([0 1 0].');   % explanatory purposes.

    for j=1:Ntheta
        PR((k-1)*Ntheta+j,:)=cos((j-1)*dtheta).*e1+sin((j-1)*dtheta).*e2;
    end

end

return;