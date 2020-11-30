function [Rest,estdx]=refine3DmatchBFGS(vol1,vol2,R1,estdx,verbose)

% Create initial guess vector
[psi,theta,phi]=rot2xyz(R1);
X0=[ psi; theta; phi; estdx(:) ];
X0=double(X0);

f = @(X)eval3Dmatchaux(X,vol1,vol2);

opt.GoalsExactAchieve=0;
opt.TolFUN=1.0e-4;
opt.TolX=1.0e-4;

if verbose
    opt.Display='iter';
else
    opt.Display='off';
end

X=fminlbfgs(f,X0,opt);

psi=X(1);
theta=X(2);
phi=X(3);
Rest=Rz(phi)*Ry(theta)*Rx(psi);
estdx=[X(4);X(5);X(6)];
