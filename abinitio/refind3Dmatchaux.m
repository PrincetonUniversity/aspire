function [Rest,estdx,optout]=refind3Dmatchaux(vol1,vol2,R1,estdx)

% Create initial guess vector
[psi,theta,phi]=rot2xyz(R1);
X0=[ psi; theta; phi; estdx(:) ];
X0=double(X0);

f = @(X)eval3Dmatchaux(X,vol1,vol2);

% options = optimoptions('fminunc','Algorithm','quasi-newton',...
%     'Display','iter-detailed','TolFun',1.0e-2,'TolX',1.0e-4);
% [X,~,~,optout]=fminunc(f,X0,options);

options = optimset('Display','iter','TolFun',1.0e-2,'TolX',1.0e-4,'MaxFunEvals',100);
[X,~,~,optout]=fminsearch(f,X0,options);


psi=X(1);
theta=X(2);
phi=X(3);
Rest=Rz(phi)*Ry(theta)*Rx(psi);
estdx=[X(4);X(5);X(6)];
