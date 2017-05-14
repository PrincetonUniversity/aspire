function [Rijs,Rijgs,Confijs] = cryo_generateRij(clmatrix,n_theta,refq)
%
%
% For each triplet of projections (i,j,k), we find the relative rotation of
% Ri and Rj using Rk. The final estimate of Rij is done using all possible
% k's that induce a feasible angle. A histogram of all feasible angles is
% constructed. Finaly, the relative-roation Rij, is constructed via the median angle
% of the angles falling withing 10 degrees away from the histogram peak
%
% input
%   clmatrix:    K-by-K-by 2, storing the indecies of commonlines c_ij and
%                c_ji, and c^g_ij and c^g_ji
%   L:           Number of Fourier lines per image.
%   refq:        (optional) are the quaternions used to compute the commonline matrix.
%
% output
%   Rijs:       9-by-(K*(K-1)/2)-by2, Rijs(:, ind,1) stores the relative
%               rotation R_i^TR_j or JR_i^TR_jJ. Rijs(:, ind,2) stores the
%               relative rotation R_i^TgR_j or JR_i^TgR_jJ where ind is the linear
%               index of (i,j), i<j, see function  uppertri_ijtoind.m
%
% Yoel Shkolnisky and Xiuyuan Cheng, July 2012.
% Modified : Gabi Pragier , November 2013

% Check the input
sz=size(clmatrix);
if numel(sz)~=3
    error('clmatrix must be a 3-d matrix');
end
if sz(1)~=sz(2)
    error('clmatrix must be a square matrix');
end
if sz(3)~=2
    error('clmatrix must be');
end

nImages=sz(1); % number of images
nPairs = nchoosek(nImages,2);
%Rij(:, (i,j)-th,1) = Ri'*Rj   up to J-conjugacy, and
%Rij(:, (i,j)-th,2) = Ri'*g*Rj up to J-conjugacy
Rijs_tmp = zeros(9, nPairs,2);
ConfijsTmp = zeros(nImages,nImages,2); % the confidence in the estimation

poolreopen;

for k1=1:nImages-1
    %fprintf('Process image %d out of %d\n',k1,K-1);
    for cl_layer=1:2
        Rijtmp=zeros(3,3,nImages);
        %for k2=k1+1:K
        parfor k2=k1+1:nImages
            [Rijtmp(:,:,k2),ConfijsTmp(k1,k2,cl_layer)] = ...
                cryo_calpha_voteIJ_wrapper(clmatrix,k1,k2,1:nImages,cl_layer,0,n_theta,refq);
        end
        for k2=k1+1:nImages
            ind =  uppertri_ijtoind(k1,k2,nImages);
            Rk = Rijtmp(:,:,k2);
            Rijs_tmp(:,ind,cl_layer) = Rk(:);
        end
    end
end

Rijs  = reshape(Rijs_tmp(:,:,1),[3,3,nPairs]);
Rijgs = reshape(Rijs_tmp(:,:,2),[3,3,nPairs]);

Confijs = mean(ConfijsTmp,3);

end


function [Rk1k2,Kfeasible] = cryo_calpha_voteIJ_wrapper(clmatrix,k1,k2,K3,cl_layer,is_perturbed,n_theta,refq)

[~, ~, cosalpha, Kfeasible]=cryo_calpha_voteIJ(clmatrix,k1,k2,K3,cl_layer,is_perturbed,n_theta,refq);

ref=0;
if exist('refq','var')
    ref=1;  % Reference quaternions are given.
end

TOL=1.0E-12; % This tolerance is relevant only if L is very large (1.0e15),
% to check for accuracy. Otherwise, ignore it and the
% resulting output.
%ZERO = 1e-16;

J=[1 0 0; 0 1 0; 0 0 -1]; % Reflection matrix


%good_k3 =cryo_voteIJ(clmatrix,L,k1,k2,1:K,refq,0);

if ~isempty(cosalpha)
    %taking the averaged cos_alpha from good k3's that is in the peak of voting
    cos_alpha_averaged = median(cosalpha);
    
    %compute R1^T*R2 (or JR1^T*R2J ) by c12, c21 and the alpha detected
    %by the voting over k3
    cl12= clmatrix([k1 k2],[k1 k2 ],cl_layer);
    
    aa=(cl12(1,2)-1)*2*pi/n_theta; %the euler angles
    bb=(cl12(2,1)-1)*2*pi/n_theta;
    
    alpha = acos(cos_alpha_averaged); %% alpha is phi21
    % GABI : note that alpha could be also -pi (acos is same for both)!!.
    %(this is the manefistation for the J-conjugacy uncertainty)
    
    % the relative rotation R12=R1^TR2 from the euler angle beta1 (on I1 from x-axis to c12),
    % phi12 (the tilting angle from I2 to I1) , beta2 (on I2 from x-axis to c21)
    % phi12 is calculated from voting (see above). R12 =
    % O_{ZXZ}(beta2, phi21, -beta1). Due to the seeting in pft,
    % beta1= pi-aa, beta2 = pi-bb
    R12 = ang2orth(pi-bb, alpha, -(pi-aa));
    
    [U, ~, V] = svd(R12);
    Rk1k2 = U*V.'; %Rk is the estimation of R1^TR2 up to J-conjugacy
    
    %             k3 = good_k3(1);
    %             clk= clmatrix([k1 k2 k3],[k1 k2 k3]);
    %             Rk_ = rotratio_eulerangle(clk,L);
    %
    %             [R12, Rk_]
    
    
    if ref
        % Compare resulting rotation computed using common
        % lines to the rotations computed using the true
        % rotations.
        R1ref=q_to_rot(refq(:,k1));
        R2ref=q_to_rot(refq(:,k2));
        inv_R1ref=R1ref.';
        inv_R2ref=R2ref.';
        
        Rk_t1 = inv_R1ref.'*inv_R2ref;
        Rk_t2 = J*inv_R1ref.'*inv_R2ref*J;
        
        err = min( norm(Rk1k2-Rk_t1, 'fro'), norm(Rk1k2-Rk_t2, 'fro'));
        
        if (err>TOL) && (n_theta>1.0e14)
            %                    fprintf('Wrong rotation: [k1=%d  k2=%d  k3=%d] err=%e  tol=%e\n',k1,k2,k3,err,TOL);
        end
    end
    
else
    % Images k1 and k2 correspond to the same viewing direction and
    % differ only by in-plane rotation. No triangle can be formed by
    % any image k3. Thus, we find the (in-plane) rotation matrix
    % between k1 and k2 by aligning the two images. The matrix is
    % [ cos(theta) sin(theta) 0 ;...
    %  -sin(theta) cos(theta) 0 ;...
    %       0           0     1 ]
    % Here we cheat and compute it using the quaternions.
    % If refq is not given, just put zero.
    Rk1k2=zeros(3,3);
    if ref
        R1ref=q_to_rot(refq(:,k1));
        R2ref=q_to_rot(refq(:,k2));
        inv_R1ref=R1ref.';
        inv_R2ref=R2ref.';
        
        % The upper-left 2x2 block of Rk is equal
        % to the upper-left 2x2 block of
        % (inv_R1ref.'*inv_R2ref+J*inv_R1ref.'*inv_R2ref*J)/2.
        Rk1k2=inv_R1ref.'*inv_R2ref;
    end
end

end