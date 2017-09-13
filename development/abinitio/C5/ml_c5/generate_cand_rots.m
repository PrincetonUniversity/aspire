function [Ris_tilde,R_theta_ijs] = generate_cand_rots(nPoints_sphere,inplane_rot_res,is_use_gt_in_cands,refq)

vis = generate_vis(nPoints_sphere,is_use_gt_in_cands,refq);

nRisTilde = size(vis,1);
Ris_tilde  = zeros(3,3,nRisTilde);
for i=1:nRisTilde
    vi = vis(i,:);
    Ris_tilde(:,:,i) = complete_3rdRow_to_rot(vi);
end


% create all in-plane angular differences using 'inplane_rot_res' as the
% discretization
theta_ij = (0:inplane_rot_res:(360-inplane_rot_res))*pi/180;
n_theta_ij = numel(theta_ij);

% TODO: this should be enforced in the interface
assert(mod(n_theta_ij,4)==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: construct all in-plane rotation matrices R(theta_ij)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cos_theta_ij = cos(theta_ij);
sin_theta_ij = sin(theta_ij);
zrs = zeros(1,n_theta_ij);
ons = ones(1,n_theta_ij);
R_theta_ijs = [cos_theta_ij  ; sin_theta_ij; zrs; ...
    -sin_theta_ij  ; cos_theta_ij; zrs; ...
    zrs;            zrs;          ons];

R_theta_ijs = reshape(R_theta_ijs,3,3,n_theta_ij);

end


function vis = generate_vis(nPoints_sphere,is_use_gt,refq)

if ~is_use_gt
    vis = randn(nPoints_sphere,3);
    for i=1:nPoints_sphere
        % normalize (TODO: what happens if norm is too small. need to discard)
        vis(i,:) = vis(i,:)./norm(vis(i,:));
    end
    
    nImages = size(refq,2);
    vis_gt = zeros(nImages,3);
    for i=1:nImages
        rot_gt = q_to_rot(refq(:,i)).';
        vis_gt(i,:) = rot_gt(3,:);
    end
    
    test_coverage(vis.',vis_gt.','red','blue');
else
    nImages = size(refq,2);
    vis_gt = zeros(nImages,3);
    for i=1:nImages
        rot_gt = q_to_rot(refq(:,i))';
        vis_gt(i,:) = rot_gt(3,:);
        
    end
    
    vis_pert = zeros(nImages,3);
    for i=1:nImages
        % add some pertrubation to each row
        noise = randn(1,3);
        noise = tand(5)*noise./norm(noise);
        vis_pert(i,:) = vis_gt(i,:) + noise;
        vis_pert(i,:) = vis_pert(i,:)./norm(vis_pert(i,:));
    end
    
    vis_bad = zeros(2*nImages,3);
    for i=1:size(vis_bad,1)
        tmp = randn(3,3); [Q,~] = qr(tmp);
        vis_bad(i,:) = Q(3,:);
    end
        
    vis = [vis_gt ; vis_pert; vis_bad];
end

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



function test_coverage(vis,vis_gt,color1,color2)

figure;
plot_view_dirs(vis,color1);
hold on;
plot_view_dirs(vis_gt,color2);
% pause;

end