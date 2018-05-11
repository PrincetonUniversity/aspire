function [Rest,estdx,optout]=optimize_orientation(proj_hat,R,refprojs_hat,Rrefs,L,estdx)
%
% Given a (Fourier transformed) projection proj_hat, estimate its rotation
% and translation parameters by finding its common lines with reference
% (Fourier transformed) projections refprojs_har whose rotations Rrefs are
% known.
%
% NOTE: proj_hat and refprojs_hat are assumed to be already normalized
% using cryo_raynormalize.

% Create initial guess vector
[psi,theta,phi]=rot2xyz(R);
X0=[ psi; theta; phi; estdx(:) ];
X0=double(X0);

% Nrefs=size(refprojs_hat,3);
% for k=1:Nrefs
%     pf=refprojs_hat(:,:,k);    
%     pf=cryo_raynormalize(pf);
%     refprojs_hat(:,:,k)=pf;
% end
% 
% proj_hat=cryo_raynormalize(proj_hat);

f = @(X)evalRmatchaux(X,proj_hat,refprojs_hat,Rrefs,L);

% options = optimoptions('fminunc','Algorithm','quasi-newton',...
%     'Display','iter-detailed','TolFun',1.0e-2,'TolX',1.0e-4);
% [X,~,~,optout]=fminunc(f,X0,options);

options = optimset('Display','off','TolFun',1.0e-4,'TolX',1.0e-4,'MaxFunEvals',200);
[X,~,~,optout]=fminsearch(f,X0,options);

psi=X(1);
theta=X(2);
phi=X(3);
Rest=Rz(phi)*Ry(theta)*Rx(psi);
estdx=[X(4);X(5)];
