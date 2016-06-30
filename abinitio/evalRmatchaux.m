function e=evalRmatchaux(X,proj_hat,refprojs_hat,Rrefs,L)

% First three parameters of X are the (psi,theta,phi) representation of the
% rotation we are optimizing over. Convert these angles back to a rotation
% matrix so we can easily extract common lijes.

psi=X(1);
theta=X(2);
phi=X(3);
R=Rz(phi)*Ry(theta)*Rx(psi);
dx=[X(4);X(5)];

% Find common lines between rotation and reference rotations.
[Mkj,Ckj,Cjk]=commonline_R_vec(Rrefs,R,L,0.999);

dtheta=2*pi/L;
rmax=size(proj_hat,1);
rk2=(0:rmax-1).';

Nrefs=size(refprojs_hat,3);
cvals=zeros(Nrefs,1);

for j=1:Nrefs
    if Mkj(j)>0
        U=proj_hat(:,Cjk(j));
        V=refprojs_hat(:,Ckj(j),j);
        V=bsxfun(@times,conj(U),V);
        
        alpha=(Cjk(j)-1)*dtheta;  % Angle of common ray in projection k.
        delta=sin(alpha)*dx(1)+cos(alpha)*dx(2);      
        shift_phases=exp(+2*pi*sqrt(-1).*rk2.*delta./(2*rmax+1)); % - sign since we shift back.
        W=real(sum(V.*shift_phases));
        cvals(j)=W;
    end
end

c=sum(cvals)./sum(cvals>0);
% Extract Fourier rays corresponding to the common lines

% Reshift Fourier rays

% Compute correlations

e=1-c;