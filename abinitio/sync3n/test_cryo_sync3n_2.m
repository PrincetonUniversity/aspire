% Test the function cryo_sync3n_estimate rotations.m
% 
% All messages should print OK, with errors of the order of machine
% precision.
%
% Yoel Shkolnisky, March 2016.

initstate;
TOL=1.0e-14;
K=50;
refq=qrand(K);
L=1E15;     % Use a large number of lines per image, so we don't have discretization errors.
cl=clmatrix_cheat_q(refq,L);

open_log(0);
rotations=cryo_sync3n_estimate_rotations(cl,L);

Rref=zeros(3,3,K);
for k=1:K
    R=q_to_rot(refq(:,k));
    Rref(:,:,k)=R.';
end

[~,mse,diff,~,~]=register_rotations(rotations,Rref);

% Check that the MSE of the recovered orientations is very small
e2=max(diff(:));

if e2>TOL
    log_message('**** Maximal rotation is error is large e2=%e\n',e2);
else
    log_message('Maximal rotation is error is OK e2=%e\n',e2);
end

if mse>TOL
    log_message('**** mse is large mse=%e\n',mse);
else
    log_message('mse OK mse=%e\n',mse);
end

