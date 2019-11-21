% Use the SDP method and experiment with fourier transformations of the
% SO(3) with more than two molecules
% 09/24/2015

%function [error] = SDP_subtomo(n, K_sub, SNR)
%% 1. Preparation of the underlying rotation matrices

% STEP 0
% Generate simulated projections of size 65x65.
% For simplicity, the projections are centered.
n = 65;
ns = 9;
K_sub = 3;
SNR = 1000; 
K = K_sub*ns;

% STEP 1
% First, rotate the volume to a random angle to represent the random
% rotation of molecule.
% Then, create the refq matrix where each column is a quaternion of the rotation,
% the first nine columns correspond to the first subtomogram, while the
% last nine columns represent the second subtomogram. The tilt degress are
% . The axis is fixed as (1, 0, 0) for each subtomogram.
% The quaternions are in the form q = cos(t/2) + (x,y,z)sin(t/2).

% rotate around the y axis such that the x-axis rotate to the -z-axis
%%%% add initial position of the second molecule Sep14, 2015 %%%%
subq = qrand(K_sub);
tiltq = zeros(4,K);   
nh = (ns-1)/2;
grid = repmat((-nh:nh)*pi/24,1,K_sub);     % tilt angles in radians, their relative
                                        % angle is 7.5 deg = 1/24*pi
tiltq(1,:) = cos(grid/2);   %%%%%% use grid/2 to fit the quaternion definition            
tiltq(2,:) = sin(grid/2);         % assume the rotation axis of tilt series 
                                % is the x-axis, the projection is in the z-axis
ref_q = zeros(4,K); % refq = b*a = q*r which is rotation b follows rotation a in absolute frame ### check
for i = 1:K
    a = subq(:,ceil(i/ns));
    b = tiltq(:,i);
    ref_q(1,i) = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4);
    ref_q(2,i) = a(1)*b(2) + a(2)*b(1) - a(3)*b(4) + a(4)*b(3);
    ref_q(3,i) = a(1)*b(3) + a(2)*b(4) + a(3)*b(1) - a(4)*b(2);
    ref_q(4,i) = a(1)*b(4) - a(2)*b(3) + a(3)*b(2) + a(4)*b(1);
end
% ref_q = subq;

%% 2. Generate images
% create projection slices with the quaternions
[~,noisy_projs,~,~]=cryo_gen_projections_vsub(n,K,SNR,ref_q);
%masked_projs=mask_fuzzy(projs,33); % Applly circular mask
masked_noisy_projs = mask_fuzzy(noisy_projs,33);

% Compute polar Fourier transform, using radial resolution n_r and angular
% resolution n_theta. 
n_theta=72;
n_r=65;
[npf,~]=cryo_pft(masked_noisy_projs,n_r,n_theta,'single');

% Find common lines from clean projections, for each corrstack (:,k,l) the
% first array indices are i*n_theta*j*n_theta/2 where the i-th slice in subtomogram A
% and the j-th j<=n_theta/2 slice in subtomogram B -> use ind2sub/sub2ind
max_shift=0;
shift_step=1;
% tic
%[clstack, ~, ~, ~] = commonlines_gaussian_vsub(npf,max_shift,shift_step);
[clstack, ~, ~, ~] = commonlines_gaussian_vsubeucl(npf,max_shift,shift_step);
%[ref_clstack,~]=clmatrix_cheat_q(ref_q,n_theta);
% toc
% disp('time of computing the common line matrix');
% disp(nnz(clstack - ref_clstack));

%% 3. build SDP coefficient matrix from common lines
dim = 2;
% incorporate the rotation of the slices in each subtomograms
Cs = zeros(dim,dim*ns);
for i = 1:ns
    R = q_to_rot(tiltq(:,i));
    Cs(:,dim*i-1:dim*i) = R(1:dim,1:dim)'; % possible
end

C = zeros(dim*K_sub);
d_ang = 2*pi/n_theta;
% iterate over subtomograms
for i = 1:K_sub-1
    for j = i+1:K_sub
        % iterate over slices
        for k = 1:ns
            for l = 1:ns
                I = (i-1)*ns+k;
                J = (j-1)*ns+l;
                n_i = d_ang*subplus(clstack(I,J)-1);
                n_j = d_ang*subplus(clstack(J,I)-1);
                C(dim*i-1:dim*i,dim*j-1:dim*j) = C(dim*i-1:dim*i,dim*j-1:dim*j) ...
                    + Cs(:,dim*l-1:dim*l)*[cos(n_j); sin(n_j)]*[cos(n_i) sin(n_i)]*Cs(:,dim*k-1:dim*k)';
            end
        end
    end
end
C = C+C';

%% SDP optimization
N = dim*K_sub;

% use cvx to optimize
cvx_begin SDP
    variable G(N,N) semidefinite; 
    maximize(trace(C*G))
    subject to
        for i = 1:K_sub
            G(dim*i-1:dim*i,dim*i-1:dim*i) == eye(dim);
        end
cvx_end

% works with Cryo-em single slices, bug in the slice indices

%% analyze
R = zeros(dim,dim*K_sub);
for i = 1:K_sub
    R_sub = q_to_rot(subq(:,i));
    R(:,dim*i-1:dim*i) = R_sub(1:dim,1:dim)';
end
G_true = R'*R;
G_trueT = eye(dim*K_sub); % column i for molecule i with others rotation
for i = 1:K_sub-1
    for j = i+1:K_sub
        G_trueT(dim*i-1:dim*i,dim*j-1:dim*j) = G_true(dim*i-1:dim*i,dim*j-1:dim*j)';
    end
end
G_trueT = G_trueT+G_trueT'-eye(dim*K_sub);
        
error = norm(G-G_trueT,'fro')/norm(G_trueT,'fro');
%end

