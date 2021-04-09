function [R_est,R_est_J] = fastAlignment3D(sym,vol1,vol2,verbose,n,N_projs,true_R,G,noise,SNR)
%% This function does the work for cryo_align_vols.    
%% input: 
% sym- the symmetry type- 'Cn'\'Dn'\'T'\'O'\'I', were n is the the symmetry
%      order.  
% vol1- 3D reference volume that vol2 should be aligned accordingly.
% vol2- 3D volume to be aligned.
% verbose- enter some number different from 0 if you wish to print log 
%       messages. (default is 0).
% n- the size of vol1 and vol2.
% N_projs- number of reference projections for the alignment. 
% true_R- the true rotation matrix between vol_2 and vol_1. 
% G_- size=(3,3,n) all n symmetry group elemnts. 
% noise- enter a number different from 0 in order to add noise to the
%       projections.
% SNR- if noise ~= 0 then enter the required SNR.  
%% output: 
% R_est- the estimated rotation between vol_2 and vol_1 without reflection.
% R_est_J- the estimated rotation between vol_2 and vol_1 with reflection.

%% 
refrot = 1;
if isempty(true_R), refrot = 0; end  

currentsilentmode = log_silent(verbose == 0);

%% Generate reference projections from vol2:
log_message('Generating %d projections from vol2',N_projs);
Rots = genRotationsGrid(75);
sz_Rots = size(Rots,3);

R_ref = Rots(:,:,randperm(sz_Rots,N_projs));       % size (3,3,N_projs)
ref_projs = cryo_project(vol2,R_ref,n);
ref_projs = permute(ref_projs,[2 1 3]);
R_ref = permute(R_ref,[2,1,3]);         % the true rotations.

%% Add noise to reference projections:
if noise ~= 0
    ref_projs = cryo_addnoise(ref_projs, SNR, 'gaussian');
    log_message('Added noise to reference projections with SNR=%5.3f.',SNR);
end
%% Align reference projections to vol1:
opt.N_projs = N_projs; opt.G = G; opt.Rots = Rots;
log_message('Aligning projections of vol2 to vol1');
if refrot == 1
    R = true_R;
    R = R.';
    R = R([2 1 3],[2 1 3]);
    true_R_tild = zeros(3,3,N_projs);
    true_R_tild_J = zeros(3,3,N_projs);
    for i = 1:N_projs
        true_R_tild(:,:,i) = (R*R_ref(:,:,i));
        J3 = diag([1 1 -1]);
        true_R_tild_J(:,:,i) = (J3*R*J3*R_ref(:,:,i));
    end   
    opt.true_Rots = true_R_tild; opt.true_Rots_J = true_R_tild_J;
    [R_tild,~] = cryo_align_projs(sym,ref_projs,vol1,verbose,opt);     % size (3,3,N_projs). 
else
    [R_tild,~] = cryo_align_projs(sym,ref_projs,vol1,verbose,opt);     % size (3,3,N_projs). 
end
%XXX The following comment wasn't not clear enough for me. XXX
%XXX Also, there is no description of O. XXX.
%% Synchronization:
% define the relation between the rotations as:
% 1. there is no reflection 
% between vol1 and vol2, thus the relation is only by rotation: Ri = O*g^{si}*Ri_tild,
% where Ri is the true rotation in the coordinates of vol2- R_ref(i), 
% and Ri_tild is the estimated rotation in the coordinates of vol 1-
% R_tild(i). and g^{si} is the symmetry group element of image i, when si
% in {0,1,...n-1}. 
% 2. there is reflection, in that case the relation can be
% written as: Ri = O*g^{si}*(J*Ri_tild*J). 
% if we define Oi = Ri*Ri_tild.' (or Oi = Ri*(J*Ri_tild*J).') than we get
% Oi.'*Oj = g^{si}.'*g^{sj}= g^{sj-si}. 
% thus, we can estimate all relative group elements and construct the 
% synchronization matrix Oij = Oi.'*Oj. then estimate the group elemnts for 
% each image with and whitout reflection, and latter choose the better option. 

%%% estimate O with and without reflection:
O_mat = zeros(3,3,N_projs);
O_mat_J = zeros(3,3,N_projs);
J3 = diag([1 1 -1]);
for i = 1:N_projs
    O_mat(:,:,i) = R_ref(:,:,i)*R_tild(:,:,i).';
    O_mat_J(:,:,i) = R_ref(:,:,i)*(J3*R_tild(:,:,i)*J3).';
end

%%% construct the synchronization matrix with and without reflection:
O_ij = zeros(3*N_projs,3*N_projs);
O_ij_J = zeros(3*N_projs,3*N_projs);
for i = 1:N_projs
    for j = i+1:N_projs
        O_ij((3*(i-1)+1):(3*i),(3*(j-1)+1):(3*j)) = O_mat(:,:,i).'*O_mat(:,:,j);
        O_ij_J((3*(i-1)+1):(3*i),(3*(j-1)+1):(3*j)) = O_mat_J(:,:,i).'*O_mat_J(:,:,j);
    end
end
% Enforce symmetry:
O_ij = O_ij+O_ij.'; O_ij = O_ij + eye(size(O_ij)); 
O_ij_J = O_ij_J+O_ij_J.'; O_ij_J = O_ij_J + eye(size(O_ij_J)); 

%%% define v = [g^{s1}.',..., g^{sN_projs}.'].' (v is of size 3*N_projx3),
%%% then Oij = v*v.', and Oij*v = N_projs*v. so v is an eigenvector of Oij. 
%%% the matrix Oij should be of rank 3. 
%%% find the top 3 eigenvectors:
% without reflection:
[U,s]=eig(O_ij); s=diag(s); [s,ii]=sort(s,'descend'); U=U(:,ii);
V = U(:,1:3);

% with reflection:
[UJ,sJ]=eig(O_ij_J); sJ=diag(sJ); [sJ,ii]=sort(sJ,'descend'); UJ=UJ(:,ii);
VJ = UJ(:,1:3);


% XXX From the printouts, both matrices are rank-3. Is this ok? XXX
%%% Print the first singular values of both synchronization matrices. Only
%%% one of them should be rank-3. In fact, at this point we can figure out if
%%% the two volumes are reflected with respect to each other or not. However,
%%% we make this decision later in the code.
log_message('First 6 singular values of synchronization matrix:');
log_message('\t without reflection (%4.2e, %4.2e, %4.2e, %4.2e, %4.2e, %4.2e)',...
        s(1),s(2),s(3),s(4),s(5),s(6));
log_message('\t with    reflection (%4.2e, %4.2e, %4.2e, %4.2e, %4.2e, %4.2e)',...
        sJ(1),sJ(2),sJ(3),sJ(4),sJ(5),sJ(6));
log_message('In the noiseless case the synchronization matrix should be rank 3.');
% XXX Why not using the rank information to detect reflection? Isn't it
% more robust than correlation? XXX


% XXX This part is not commented. What is G and what is the purpose of this
% code? XXX
%%% estimating G:
G = zeros(3,3,N_projs);
G_J = zeros(3,3,N_projs);
for i = 1:N_projs
    B = V((3*(i-1)+1):(3*i),:);
    [u_tmp,~,v_tmp] = svd(B);
    B_round = det(u_tmp*v_tmp.').*u_tmp*v_tmp.';
    G(:,:,i) = B_round.';
    
    BJ = VJ((3*(i-1)+1):(3*i),:);
    [uJ_tmp,~,vJ_tmp] = svd(BJ);
    BJ_round = det(uJ_tmp*vJ_tmp.').*uJ_tmp*vJ_tmp.';
    G_J(:,:,i) = BJ_round.';
end

%%% set the global rotation to be an element from the symmetry group: 
O1 = G(:,:,1).';
O1_J = G_J(:,:,1).';
for i=1:N_projs
    G(:,:,i) = O1*G(:,:,i);
    G_J(:,:,i) = O1_J*G_J(:,:,i);
end

%% Estimating the rotation:
%%% Estimate the two candidate orthogonal transformations.
%%% Only one of them would be orthogonal, depending if we have or not
%%% reflection between the volumes.
O_mat1 = zeros(3,3,N_projs);
O_mat_J1 = zeros(3,3,N_projs);
for i = 1:N_projs
    O_mat1(:,:,i)=O_mat(:,:,i)*(G(:,:,i)).';
    O_mat_J1(:,:,i)=O_mat_J(:,:,i)*(G_J(:,:,i)).';
end

O = mean(O_mat1,3);
O_J = mean(O_mat_J1,3);

R = O;
[U,~,V] = svd(R); % Project R to the nearest rotation.
R_est = U*V.';
assert(det(R_est)>0);
R_est = R_est([2 1 3],[2 1 3]);
R_est = R_est.';

R_J = O_J;
[U,~,V] = svd(R_J); % Project R to the nearest rotation.
R_est_J = U*V.';
assert(det(R_est_J)>0);
R_est_J = R_est_J([2 1 3],[2 1 3]);
R_est_J = R_est_J.';

log_silent(currentsilentmode);
end
