nGT_rots = 20;
snr = 10000000;
max_shift = 0;
shift_step =1;
n_theta = 360;

nPoints = 400;
inplane_rot_res = 1;
is_test_coverage = true;

[~,refq] = generate_c4_images(nGT_rots,snr,65,'GAUSSIAN',max_shift,shift_step);

rots_gt = zeros(3,3,nGT_rots);
for i=1:nGT_rots
    rots_gt(:,:,i) = q_to_rot(refq(:,i)).';
end

% sample points on sphere and construct Ri_tilde
vis = randn(nPoints,3);
Ris_tilde = zeros(3,3,nPoints);
for i =1:nPoints
    vi = vis(i,:);
    vi = vi./norm(vi);
    Ris_tilde(:,:,i) = complete_3rdRow_to_rot(vi);
end

if is_test_coverage
    vis = squeeze(Ris_tilde(3,:,:));
    vis_gt = squeeze(rots_gt(3,:,:));
    test_coverage(vis,vis_gt,'red','blue');
end

% create all in-plane angular differences using 'inplane_rot_res' as the
% discretization
theta_ij = (0:inplane_rot_res:(360-inplane_rot_res))*pi/180;
n_theta_ij = numel(theta_ij);

cos_theta_ij = cos(theta_ij);
sin_theta_ij = sin(theta_ij);
zrs = zeros(1,n_theta_ij);
ons = ones(1,n_theta_ij);
R_theta_ij = [cos_theta_ij  ; sin_theta_ij; zrs; ...
    -sin_theta_ij  ; cos_theta_ij; zrs; ...
    zrs;            zrs;          ons];

R_theta_ij = reshape(R_theta_ij,3,3,n_theta_ij);


rots_cands = zeros(3,3,nPoints,n_theta_ij);
for i=1:nPoints
    for j=1:n_theta_ij
        rots_cands(:,:,i,j) = R_theta_ij(:,:,j)*Ris_tilde(:,:,i);
    end
end
rots_cands = reshape(rots_cands,[3,3,nPoints*n_theta_ij]);
nCands = size(rots_cands,3);

% find for each of the GT rots which is the closest candidate rotation
min_dists    = zeros(1,nGT_rots);
closest_cand = zeros(1,nGT_rots);
est_rots = zeros(3,3,nGT_rots);
msg = [];
for i=1:nGT_rots
    t1 = clock;
    dist = zeros(1,nCands);
    for j=1:nCands
        dist(j) = norm(rots_gt(:,:,i)-rots_cands(:,:,j),'fro');
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
    msg = sprintf('i=%3d/%3d j=%3d/%3d t=%7.5f',i,nGT_rots,j,nCands,t);
    fprintf('%s',msg);
    %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


log_message('median view direction angular error=%.2f',median(acosd(1-min_dists)));
if is_test_coverage
    vis = squeeze(est_rots(:,3,:));
    vis_gt = squeeze(rots_gt(:,3,:));
    test_coverage(vis,vis_gt,'red','blue');
end




nImages = size(refq,2);
assert(size(est_rots,3) == size(refq,2));
d_theta = 2*pi/n_theta;
err_in_degrees = zeros(1,n_theta*nImages);
for k = 1:nImages
    rays_gt = zeros(n_theta,3);
    Rk_gt   = q_to_rot(refq(:,k))';
    for j = 1:n_theta
        rays_gt(j,:) = cos((j-1)*d_theta).*Rk_gt(:,1) ...
            +  ...
            sin((j-1)*d_theta).*Rk_gt(:,2);
    end
    rays = zeros(n_theta,3);
    Rk   = est_rots(:,:,k);
    for j = 1:n_theta
        rays(j,:) = cos((j-1)*d_theta).*Rk(:,1) ...
            +  ...
            sin((j-1)*d_theta).*Rk(:,2);
    end
    
    cos_angle   = rays(:,1).*rays_gt(:,1) + rays(:,2).*rays_gt(:,2) + rays(:,3).*rays_gt(:,3);
    err         = acos(cos_angle);
    err_in_degrees((k-1)*n_theta+1:k*n_theta) = err * 180/pi;
end

log_message('\nmedian error in degrees: %e',median(err_in_degrees));
log_message('mean error in degrees: %e',mean(err_in_degrees));
log_message('std error in degrees: %e',std(err_in_degrees));
log_message('min error in degrees: %e',min(err_in_degrees));
log_message('max error in degrees: %e',max(err_in_degrees));
