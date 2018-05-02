function [mse1, mse2] = stam

nGT_rots = 30;
nPoints_sphere  = 500;
inplane_rot_res = 1;
snr = 100000000000000;
max_shift = 0;
shift_step = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate fround truth rotation matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,refq] = generate_c4_images(nGT_rots,snr,65,'GAUSSIAN',max_shift,shift_step);

rots_gt = zeros(3,3,nGT_rots);
for i=1:nGT_rots
    rots_gt(:,:,i) = q_to_rot(refq(:,i)).';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample sphere points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vis = randn(nPoints_sphere,3);
for i=1:nPoints_sphere
    % normalize (TODO: what happens if norm is too small. need to discard)
    vis(i,:) = vis(i,:)./norm(vis(i,:));
end


theta_ij = (0:inplane_rot_res:(360-inplane_rot_res))*pi/180;
n_theta_ij = numel(theta_ij);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct all in-plane rotation matrices R(theta_ij)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cos_theta_ij = cos(theta_ij);
sin_theta_ij = sin(theta_ij);
zrs = zeros(1,n_theta_ij);
ons = ones(1,n_theta_ij);
R_theta_ijs = [cos_theta_ij  ; sin_theta_ij; zrs; ...
              -sin_theta_ij  ; cos_theta_ij; zrs; ...
                    zrs;            zrs;     ons];

R_theta_ijs = reshape(R_theta_ijs,3,3,n_theta_ij);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method 1: as in code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ris_tilde_1  = zeros(3,3,nPoints_sphere);
for i=1:nPoints_sphere
    vi = vis(i,:);
    Ris_tilde_1(:,:,i) = complete_3rdRow_to_rot(vi);
end

R_theta_ijs_1 = R_theta_ijs;
R_theta_ijs_1(:,:,ceil(3/4*end+1):end) = [];
rots_cands_1 = form_Ris(Ris_tilde_1,R_theta_ijs_1,1);

mse1 = match(rots_gt,rots_cands_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method 2: as per geometric intuition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_theta_ijs_2 = R_theta_ijs;
beam_dirs = vis';
beam_dirs_azimut = atan2(beam_dirs(2,:),beam_dirs(1,:));
inds_first_qtr = find(beam_dirs_azimut > 0 & beam_dirs_azimut < pi/2);

beam_dirs = beam_dirs(:,inds_first_qtr);
nbeam_dirs = size(beam_dirs,2);
Ris_tilde_2  = zeros(3,3,nbeam_dirs);
for i=1:nbeam_dirs
    beam_dir = beam_dirs(:,i);
    Ris_tilde_2(:,:,i) = transpose(complete_3rdRow_to_rot(beam_dir'));
end

rots_cands_2 = form_Ris(Ris_tilde_2,R_theta_ijs_2,2);
mse2 = match(rots_gt,rots_cands_2);

end


function R = complete_3rdRow_to_rot(r3)
%
% Constructs a rotation matrix whose third row is equal to a given row vector
% 
% Input parameters:
%   r3         A 1x3 vector of norm 1
% Output parameters:
%   R          A rotation matrix whose third row is equal to r3

assert(abs(norm(r3)-1)<1e-5); 
% handle the case that the third row coincides with the z-axis
if norm(r3-[0,0,1]) < 1e-5
    R = eye(3);
    return;
end

% tmp is non-zero since r3 does not coincide with the z-axis
tmp = sqrt(r3(1)^2 + r3(2)^2);
% contruct an orthogonal row vector of norm 1
r1 = [r3(2)/tmp -r3(1)/tmp 0];
% construct r2 so that r3=r1xr2
r2 = [r3(1)*r3(3)/tmp r3(2)*r3(3)/tmp -tmp];

R = [r1; r2; r3];

end


function Ris = form_Ris(Ris_tilde,R_theta_ijs,method_type)

nPoints_sphere = size(Ris_tilde,3);
n_theta_ij = size(R_theta_ijs,3);

Ris = zeros(3,3,nPoints_sphere*n_theta_ij);
counter = 0;
for i=1:nPoints_sphere
    for j=1:n_theta_ij
        counter = counter + 1;
        
        Ri_tilde = Ris_tilde(:,:,i);
        R_theta_ij = R_theta_ijs(:,:,j);
        if method_type == 1
            Ris(:,:,counter) = R_theta_ij*Ri_tilde;
        else
            Ris(:,:,counter) = Ri_tilde*R_theta_ij;
        end
    end
end

end

function mse = match(rots_gt,rots_cands)

g = [cosd(90) -sind(90) 0; 
	 sind(90)  cosd(90) 0; 
	 0 				 0  1]; % a rotation of 90 degrees about the z-axis

nGT_rots    = size(rots_gt,3);
nCands_rots = size(rots_cands,3);
% find for each of the GT rots which is the closest candidate rotation
min_dists    = zeros(1,nGT_rots);
closest_cand = zeros(1,nGT_rots);
est_rots = zeros(3,3,nGT_rots);
msg = [];
for i=1:nGT_rots
    t1 = clock;
    dist = zeros(1,nCands_rots);
    for j=1:nCands_rots
        dist_tmp = zeros(1,4);
        for s=0:3
            dist_tmp(s+1) = norm(rots_gt(:,:,i)-g^s*rots_cands(:,:,j),'fro');
        end
        dist(j) = min(dist_tmp);
%         dist(j) = 1-dot(rots_gt(:,3,i),rots_cands(:,3,j)); % take 1-dot so that min is optimal
    end
    [YY,II] = min(dist);
    min_dists(i) = YY;
    closest_cand(i) = II;
    est_rots(:,:,i) = rots_cands(:,:,II);
    %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t2 = clock;
    t = etime(t2,t1);
    bs = char(repmat(8,1,numel(msg)));
    fprintf('%s',bs);
    msg = sprintf('i=%3d/%3d t=%7.5f',i,nGT_rots,t);
    fprintf('%s',msg);
    %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
mse = mean(min_dists);
log_message('mse=%.2f',mse);

end
