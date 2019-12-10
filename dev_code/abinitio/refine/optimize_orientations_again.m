function [estR2,estdx2,optout]=optimize_orientations_again(projs,estR,estdx,L,refR,refdx)
%
% Given a (Fourier transformed) projection proj_hat, estimate its rotation
% and translation parameters by finding its common lines with reference
% (Fourier transformed) projections refprojs_har whose rotations Rrefs are
% known.

 % Check that all variables are of the correct dimension
 % XXXX
 
Nprojs=size(projs,3);
% Create initial guess vector;
 
eulerangs=zeros(Nprojs,3); % Each row is (psi,theta,phi) of one of the rotations
for k=1:Nprojs
    [psi,theta,phi]=rot2xyz(estR(:,:,k));
    eulerangs(k,:)=[psi, theta, phi];
end
 
X0=[ eulerangs(:); estdx(:) ];
X0=double(X0);
 
 n_r=ceil(size(projs,1)/2);
 projs_hat=cryo_pft(projs,n_r,L,'single');
 projs_hat=cryo_raynormalize(projs_hat);

f = @(X)evalmatchstackaux(X,projs_hat,L);

% options = optimoptions('fminunc','Algorithm','quasi-newton',...
%     'Display','iter-detailed','TolFun',1.0e-2,'TolX',1.0e-4);
% [X,~,~,optout]=fminunc(f,X0,options);

tol=2*pi/L;
options = optimset('Display','iter','TolFun',tol,'TolX',tol);
[X,~,~,optout]=fminsearch(f,X0,options);


% Unpack rotations
estR2=zeros(3,3,Nprojs); % Each row is (psi,theta,phi) of one of the rotations
eulerangs=reshape(X(1:3*Nprojs),Nprojs,3);
for k=1:Nprojs
    psi=eulerangs(k,1);
    theta=eulerangs(k,2);
    phi=eulerangs(k,3);
    estR2(:,:,k)=Rz(phi)*Ry(theta)*Rx(psi);
end

estdx2=zeros(2,Nprojs);
for k=1:Nprojs
    estdx2(1,k)=X(3*Nprojs+2*(k-1)+1);
    estdx2(2,k)=X(3*Nprojs+2*(k-1)+2);
end
