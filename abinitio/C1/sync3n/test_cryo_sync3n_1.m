% Test the 3Nx3N synchronization functions.
%
% The script tests the individual functions used by the 3Nx3N
% synchronization algorithm.
% 
% All messages should print OK, with errors of the order of machine
% precision.
%
% For more information see ``A graph partitioning approach to simultaneous
% angular reconstitution'', G. Pragier, I. Greenberg, X. Cheng, and Y.
% Shkolnisky, IEEE Transaction on Computational Imaging. 
%
% Yoel Shkolnisky, March 2016.


initstate;
TOL=1.0e-14;
K=50;
rots_ref = rand_rots(K);
L=1E15;     % Use a large number of lines per image, so we don't have discretization errors.
cl=clmatrix_cheat(rots_ref,L);

open_log(0);
[Rijs, ~, ~, ~] = cryo_sync3n_estimate_all_Rijs(cl, L);
        
n_eigs=3; % Number of eigenvalues of the J-sync matrix to compute.
scores_as_entries = 0; % Construct the J-sync matrix using only +1 and -1.
verbose=1;
% Estimate the handedeness of each relative rotation (up a global common
% unknown hand).
[J_sync,J_significance,eigenvalues,itr,dd] =...
    cryo_sync3n_Jsync_power_method(Rijs,n_eigs,scores_as_entries,verbose );

% Flip the handedness of the estimates Rij according to hand estimated in
% J_sync. Estimates Rij from which J_sync==-1 are J-conjugated.
Rijs = cryo_sync3n_flip_handedness(J_sync, Rijs);

% Build 3KX3K synchronization matrix 
S = cryo_sync3n_syncmatrix(Rijs);

s=eig(S);
s=sort(s);

% There should be only three large eigenvalue, and all the others very
% close to zero.
e1=norm(s(1:end-3)/norm(S)); 
if e1>TOL
    fprintf('**** Rank of S is larger than three e1=%e\n',e1);
else
    fprintf('Rank test OK. e1=%e\n',e1);
end

[Rs, ~, ~] = cryo_sync3n_S_to_rot (S, 10, 1.0e-13);

Rref=zeros(3,3,K);
for k=1:K
    R=rots_ref(:,:,k);
    Rref(:,:,k)=R.';
end

[~,mse,diff,~,~]=register_rotations(Rs,Rref);

% Check that the MSE of the recovered orientations is very small
e2=max(diff(:));

if e2>TOL
    log_message('**** Maximal rotation is error is large e2=%e',e2);
else
    log_message('Maximal rotation is error is OK e2=%e',e2);
end

if mse>TOL
    log_message('**** mse is large mse=%e',mse);
else
    log_message('mse OK mse=%e',mse);
end

