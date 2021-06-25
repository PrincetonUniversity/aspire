function [rots_alligned,err_in_degrees,mse] = analyze_results_TO(rots,symmetry,n_theta,gt_R)

[mse, rots_alligned, sign_g_Ri] = check_rotations_error(rots,gt_R,symmetry);
fprintf('MSE of rotations estimate: %f\n',mse);

[view_dir_err_in_degrees,polar_ang_gt,stats_g] = check_viewdir_degrees_error(rots_alligned,symmetry,gt_R);
fprintf('view_dir_err_in_degrees: %f\n',median(view_dir_err_in_degrees));

[err_in_degrees,mean_err_per_rot] = check_degrees_error(rots_alligned,symmetry,n_theta,gt_R);
fprintf('median error in degrees: %f\n',median(err_in_degrees));
fprintf('mean error in degrees: %f\n',mean(err_in_degrees));
fprintf('std error in degrees: %f\n',std(err_in_degrees));
fprintf('min error in degrees: %f\n',min(err_in_degrees));
fprintf('max error in degrees: %f\n',max(err_in_degrees));
end

function [mse, rot_alligned, sign_g_Ri] = check_rotations_error(est_R,gt_R,symmetry)
% Our estimate for each rotation matrix Ri may be g_iR_i, where g_i is arbitrary,
% independently of other rotation matrices. 
% As such, for error analysis we perform a g-synchronization.

sign_g_Ri = g_sync(est_R,gt_R,symmetry);

nImages = size(est_R,3);
% Register estimated rotations to true ones, and compute the difference
% between the two.
rot  = zeros(3*nImages,3); % The K estimated rotation^T matrices stacked as a matrix with 3K rows.
rot1 = zeros(3*nImages,3); % True true K rotation^T matrices stacked as a matrix with 3K rows.
rot2 = zeros(3*nImages,3); % True true K rotation^T matrices stacked as a matrix with 3K rows.

J = diag([1 1 -1]);        % Reflection matrix
gR = cryo_TO_group_elements(symmetry);

for i=1:nImages
    est_Ri = est_R(:,:,i);
    rot(3*(i-1)+1:3*i,:) = est_Ri;
    gt_Ri_gi = gt_R(:,:,i)*gR(:,:,sign_g_Ri(i)); %Rt
    rot1(3*(i-1)+1:3*i,:) = gt_Ri_gi;
    rot2(3*(i-1)+1:3*i,:) = J*gt_Ri_gi*J;
end

% Compute the two possible orthogonal matrices which register the
% estimated rotations to the true ones.
O1 = rot1.'*rot./nImages;
O2 = rot2.'*rot./nImages;

% We are registering one set of rotations (the estimated ones) to
% another set of rotations (the true ones). Thus, the transformation
% matrix between the two sets of rotations should be orthogonal. This
% matrix is either O1 if we recover the non-reflected solution, or O2,
% if we got the reflected one. In any case, one of them should be
% orthogonal.
err1 = norm(O1*O1.'-eye(3));
err2 = norm(O2*O2.'-eye(3));

errd1 = abs(abs(det(O1))-1); %allow O1 to be a reflection, i.e. in O(3)
errd2 = abs(abs(det(O2))-1); %allow O1 to be a reflection, i.e. in O(3)

% In any case, enforce the registering matrix O to be a rotation.
if err1 < err2
    [U, ~, V] = svd(O1); % Use O1 as the registering matrix
    flag=1;
else
    [U, ~, V] = svd(O2); % Use O2 as the registering matrix
    flag = 2;
end

O = V*U.';
O = O.';
if flag == 2
    O=J*O*J;
end

rot_alligned = zeros(3,3,nImages);

% Plot estimation errors
diff = zeros(nImages,1);
for i=1:nImages
    gt_Ri = gt_R(:,:,i);
    
    est_Ri= est_R(:,:,i);
    if flag ==2
        est_Ri = J*est_Ri*J;
    end
    est_Ri= est_Ri*O.'*gR(:,:,sign_g_Ri(i)).';
    rot_alligned(:,:,i) = est_Ri;
    
    diff(i) = norm(est_Ri - gt_Ri,'fro');
end

mse= sum(diff.^2)/nImages;

end

function [view_dir_err_in_degrees,polar_ang_gt,stats_g] = check_viewdir_degrees_error(rots,symmetry,gt_R)

[gR, n_gR] = cryo_TO_group_elements(symmetry);
 
nImages = size(gt_R,3);
assert(size(rots,3) == nImages);

view_dir_err_in_degrees = zeros(1,nImages);
stats_g = zeros(1,nImages);
polar_ang_gt = zeros(1,nImages);
for i = 1:nImages
    
    Ri_gt  = gt_R(:,:,i);
    qi_gt = Ri_gt(:,3);
    
    polar_ang_gt(i) = acosd(qi_gt(3));
    
    qi = rots(:,3,i);
    qis = zeros(3,n_gR);
    for s=1:n_gR
       qis(:,s) = gR(:,:,s)*qi; 
    end
    
    
    ang_diffs_g = acosd(dot(repmat(qi_gt,1,n_gR),qis));
    [YY,II] = min(ang_diffs_g);
    view_dir_err_in_degrees(i) = YY;
    stats_g(i) = II-1;
end
 
end

function [err_in_degrees,mean_err_per_rot] = check_degrees_error(rots,symmetry,n_theta,gt_R)

nImages = size(gt_R,3);
assert(size(rots,3) == nImages);
d_theta = 2*pi/n_theta;
[gR, n_gR] = cryo_TO_group_elements(symmetry);

err_in_degrees = zeros(1,n_theta*nImages);
mean_err_per_rot = zeros(1,nImages);
for k = 1:nImages
    rays_gt = zeros(n_theta,3);
    Rk_gt  = gt_R(:,:,k);
    for j = 1:n_theta
        rays_gt(j,:) = cos((j-1)*d_theta).*Rk_gt(:,1) ...
            +  ...
            sin((j-1)*d_theta).*Rk_gt(:,2);
    end
    
    c_err_in_degrees = cell(1,n_gR);
    for s = 1:n_gR
        rays = zeros(n_theta,3);
        Rk   = gR(:,:,s)*rots(:,:,k);
        for j = 1:n_theta
            rays(j,:) = cos((j-1)*d_theta).*Rk(:,1) ...
                +  ...
                sin((j-1)*d_theta).*Rk(:,2);
        end
        
        cos_angle   = rays(:,1).*rays_gt(:,1) + rays(:,2).*rays_gt(:,2) + rays(:,3).*rays_gt(:,3);
        err         = acos(cos_angle);
        c_err_in_degrees{s} = err * 180/pi;
    end
    meanErrs   = cellfun(@mean, c_err_in_degrees);
    [~, opt_s] = min(meanErrs);
    err_in_degrees((k-1)*n_theta+1:k*n_theta) = c_err_in_degrees{opt_s};
    mean_err_per_rot(k) = mean(c_err_in_degrees{opt_s});
end

[~,II] = sort(mean_err_per_rot,'descend');
n = numel(II);
pol_ang = zeros(1,n);
for i=1:n
    Ri_gt = gt_R(:,:,i);
    pol_ang(i) = acosd(Ri_gt(3,3));
end

end


function sign_g_Ri = g_sync(est_R,gt_R,symmetry)

% The in-plain rotation stage retrieves all rotation matrices. However, every calculated rotation might 
% be a the rotated version of g^s, s=1,...,n from the ground truth. 
% Moreover, the rotation to be applied might be differernt for every image.
% We therefore synchronize all rotation so that only a possibly single global
% rotation should be applied to all rotation. Since the molecule is symmetric this doesn't 
% matter as far as the reconstruction is concerned, but it is important if 
% wish to compare to ground truth

n_images = size(est_R,3);
assert(n_images == size(gt_R,3));

[gR, n_gR] = cryo_TO_group_elements(symmetry);

J = diag([1 1 -1]); % Reflection matrix

A_g = zeros(3*n_images,3*n_images);

for i=1:n_images-1 
    for j=i+1:n_images
        est_Ri = est_R(:,:,i);
        est_Rj = est_R(:,:,j);
        est_Rij = est_Ri*est_Rj.';        
        
        gt_Ri = gt_R(:,:,i);
        gt_Rj = gt_R(:,:,j);

        diffs = zeros(1,n_gR);
        for k=1:n_gR
            gt_Rij = gt_Ri*gR(:,:,k)*gt_Rj.';
            diffs(k) = min(norm(est_Rij-gt_Rij,'fro'),norm(est_Rij-J*gt_Rij*J,'fro'));
        end
        [~,ind] = min(diffs);

        A_g(3*i-2:3*i,3*j-2:3*j) = gR(:,:,ind);
    end
end

A_g = A_g + A_g';
A_g = A_g + eye(size(A_g)); 

[v,d] = eigs(A_g, 6, 'lm');
[evals,ind] = sort(diag(d), 'descend');
evect(:,1) = v(:,ind(1));
evect(:,2) = v(:,ind(2));
evect(:,3) = v(:,ind(3));

fprintf('Outer g sync first 6 eigenvalues=[%.2f %.2f %.2f %.2f %.2f %.2f]\n',...
    evals(1),evals(2),evals(3),evals(4),evals(5),evals(6));

sign_g_Ri = zeros(n_images,1);
for l=1:3
    for i=1:n_images
        vi = evect(3*i-2:3*i,l);
        vi = vi/norm(vi);           % each row should be of rank-1
        evect(3*i-2:3*i,l) = vi.';
    end
end

g_Ri = [reshape(evect(:,1),3,1,n_images) reshape(evect(:,2),3,1,n_images) reshape(evect(:,3),3,1,n_images)];

OR1 = g_Ri(:,:,1).';
for i=1:n_images
    g_Ri(:,:,i) = g_Ri(:,:,i)*OR1;
end

for i=1:n_images
    rotsDist = zeros(1,n_gR);
    for k=1:n_gR
        rotsDist(k) = norm(g_Ri(:,:,i)-gR(:,:,k),'fro');
    end
    [~,min_ind] = min(rotsDist);
    sign_g_Ri(i) = min_ind;
end 
end
