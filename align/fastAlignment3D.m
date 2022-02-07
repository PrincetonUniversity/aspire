function [R_est,R_est_J] = fastAlignment3D(sym,vol1,vol2,verbose,n,N_projs,true_R,refrot,G_group)
%% This function does the work for cryo_align_vols.    
%% input: 
% sym- the symmetry type- 'Cn'\'Dn'\'T'\'O'\'I', where n is the the 
%      symmetry order (for example: 'C2').
% vol1- 3D reference volume that vol2 should be aligned accordingly.
% vol2- 3D volume to be aligned.
% verbose- Set verbose to nonzero for verbose printouts (default is zero).
% n- the size of vol1 and vol2.
% N_projs- number of reference projections for the alignment. 
% true_R- the true rotation matrix between vol2 and vol1. 
% refrot- indicator for true_R. If true_R exist then refrot=1, else
%         refrot=0.
% G_group- size=(3,3,n) all n symmetry group elemnts.  
%% output: 
% R_est- the estimated rotation between vol_2 and vol_1 without reflection.
% R_est_J- the estimated rotation between vol_2 and vol_1 with reflection.
%% 
currentsilentmode = log_silent(verbose == 0);
%% Generate reference projections from vol2:
log_message('Generating %d projections from vol2',N_projs);
Rots = genRotationsGrid(75);
sz_Rots = size(Rots,3);
R_ref = Rots(:,:,randperm(sz_Rots,N_projs));       % size (3,3,N_projs)
ref_projs = cryo_project(vol2,R_ref,n);
ref_projs = permute(ref_projs,[2 1 3]);
R_ref = permute(R_ref,[2,1,3]);                    % the true rotations.
%% Align reference projections to vol1:
%opt.N_ref = N_projs;
opt.G = G_group; opt.Rots = Rots; opt.sym = sym;
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
    [R_tild,~] = cryo_align_projs(ref_projs,vol1,verbose,opt); % size (3,3,N_projs). 
else
    [R_tild,~] = cryo_align_projs(ref_projs,vol1,verbose,opt); % size (3,3,N_projs). 
end
%% Synchronization:
% A synchronization algorithm is used In order to revel the symmetry 
% elements of the reference projections. The relation between the volumes 
% is V2(r)=V1(Or). Denote the rotation between the volumes as X. 
% 1. In the case there is no reflection between the volumes, the rotation
%    Ri_tilde estimates giORi, therefore the approximation is O =
%    gi.'Ri_tildeRi.', where g_i is the symmetry group element of reference 
%    image i. If we define Xi=Ri*Ri_tilde.' then we get Xi.'*Xj=g_i*g_j.'. 
% 2. In the case there is a reflection between the volumes, the rotation 
%    Ri_tilde estimates qiJXRiJ, where O=JX. We have that qiJ=Jqi_tilde, 
%    therefore the approximation is X=qi_tilde.'JRi_tildeJRi.', where
%    qi_tilde is a symmetry element in the symmetry group of
%    V1_tilde(r)=V1(Jr). If we define  Xi=Ri*(J*Ri_tild*J).', then we also 
%    get Xi.'*Xj=qi_tilde*qj_tilde.'.
% Therefore, we can construct the synchronization matrix Xij=Xi.'*Xj for 
% both cases. Then, estimate the group elemnts for each image with and 
% whitout reflection, and latter choose the option that best describes the 
% relation between the two volumes. 
% estimate X with and without reflection:
X_mat = zeros(3,3,N_projs);
X_mat_J = zeros(3,3,N_projs);
J3 = diag([1 1 -1]);
for i = 1:N_projs
    X_mat(:,:,i) = R_ref(:,:,i)*R_tild(:,:,i).';
    X_mat_J(:,:,i) = R_ref(:,:,i)*(J3*R_tild(:,:,i)*J3).';
end
% construct the synchronization matrix with and without reflection:
X_ij = zeros(3*N_projs,3*N_projs);
X_ij_J = zeros(3*N_projs,3*N_projs);
for i = 1:N_projs
    for j = i+1:N_projs
        X_ij((3*(i-1)+1):(3*i),(3*(j-1)+1):(3*j)) = X_mat(:,:,i).'*X_mat(:,:,j);
        X_ij_J((3*(i-1)+1):(3*i),(3*(j-1)+1):(3*j)) = X_mat_J(:,:,i).'*X_mat_J(:,:,j);
    end
end
% Enforce symmetry:
X_ij = X_ij+X_ij.'; X_ij = X_ij + eye(size(X_ij)); 
X_ij_J = X_ij_J+X_ij_J.'; X_ij_J = X_ij_J + eye(size(X_ij_J)); 
% Define v=[g_1,..., g_N_projs].' (v is of size 3*N_projx3), then 
% X_ij=v*v.', and Xij*v=N_projs*v. Thus, v is an eigenvector of Xij. The
% matrix Xij should be of rank 3. find the top 3 eigenvectors:
% without reflection:
[U,s]=eig(X_ij); s=diag(s); [~,ii]=sort(s,'descend'); U=U(:,ii);
V = U(:,1:3);
% with reflection:
[UJ,sJ]=eig(X_ij_J); sJ=diag(sJ); [~,ii]=sort(sJ,'descend'); UJ=UJ(:,ii);
VJ = UJ(:,1:3);
% estimating G:
% Estimate the group elemnts for each reference image. G denotes the 
% estimated group without reflection, and G_J with reflection. This 
% estimation is being done from the eigenvector v by using a rounding 
% algorithm over SO(3) for each 3x3 block of v.
G = zeros(3,3,N_projs); G_J = zeros(3,3,N_projs);
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
% The global rotation from the synchronization can be any rotation matrix 
% from SO(3). So, in order to get the estimated symmetry elements to be 
% from the symmetry group we set the global rotation to be also an element 
% from the symmetry group: 
O1 = G(:,:,1).'; O1_J = G_J(:,:,1).';
for i=1:N_projs
    % setting the global rotation to be g_1.':
    G(:,:,i) = O1*G(:,:,i); G_J(:,:,i) = O1_J*G_J(:,:,i);
end
%% Estimating the rotations:
% Estimate the two candidate orthogonal transformations.
for i = 1:N_projs
    X_mat(:,:,i) = X_mat(:,:,i)*(G(:,:,i)).';
    X_mat_J(:,:,i) = X_mat_J(:,:,i)*(G_J(:,:,i)).';
end
X = mean(X_mat,3); X_J = mean(X_mat_J,3);
% Without reflection:
R = X;
[U,~,V] = svd(R); % Project R to the nearest rotation.
R_est = U*V.';
assert(det(R_est)>0);
R_est = R_est([2 1 3],[2 1 3]);
R_est = R_est.';
% with reflection:
R_J = X_J;
[U,~,V] = svd(R_J); % Project R to the nearest rotation.
R_est_J = U*V.';
assert(det(R_est_J)>0);
R_est_J = R_est_J([2 1 3],[2 1 3]);
R_est_J = R_est_J.';
log_silent(currentsilentmode);
end