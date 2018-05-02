function [rot_alligned,err_in_degrees,mse] = analyze_results_ml(rots,n_theta,refq)

[mse, rot_alligned, sign_g_Ri] = check_rotations_error(rots,refq);
log_message('MSE of rotations estimate: %e',mse);

[view_dir_err_in_degrees,polar_ang_gt,stats_g] = check_viewdir_degrees_error(rot_alligned,refq);
log_message('view_dir_err_in_degrees: %e',median(view_dir_err_in_degrees));



[err_in_degrees,mean_err_per_rot] = check_degrees_error(rot_alligned,n_theta,refq);
log_message('median error in degrees: %e',median(err_in_degrees));
log_message('mean error in degrees: %e',mean(err_in_degrees));
log_message('std error in degrees: %e',std(err_in_degrees));
log_message('min error in degrees: %e',min(err_in_degrees));
log_message('max error in degrees: %e',max(err_in_degrees));

end


function [mse, rot_alligned, sign_g_Ri] = check_rotations_error(rotations,refq)
% Our estimate for each rotation matrix Ri may be R_i, gRi, g^{2}Ri or
% g^{3}Ri independently of other rotation matrices. 
% As such, for error analysis we perform a g-synchronization.

sign_g_Ri = g_sync(rotations,refq);

nImages   = size(rotations,3);
% Register estimated rotations to true ones, and compute the difference
% between the two.
rot  = zeros(3*nImages,3);  % The K estimated rotation^T matrices stacked as a matrix with 3K rows.
rot1 = zeros(3*nImages,3); % True true K rotation^T matrices stacked as a matrix with 3K rows.
rot2 = zeros(3*nImages,3); % True true K rotation^T matrices stacked as a matrix with 3K rows.

J = diag([1 1 -1]); % Reflection matrix
g = [cosd(90) -sind(90) 0; 
	 sind(90)  cosd(90) 0; 
	 0 				 0  1]; % a rotation of 72 degrees about the z-axis
for k=1:nImages
    R = rotations(:,:,k);
    rot(3*(k-1)+1:3*k,:) = R.';
    inv_Rref = g^sign_g_Ri(k)*(q_to_rot(refq(:,k)).'); %Rt
    rot1(3*(k-1)+1:3*k,:) = inv_Rref.';
    rot2(3*(k-1)+1:3*k,:) = J*inv_Rref.'*J;
end

% Compute the two possible orthogonal matrices which register the
% estimated rotations to the true ones.
O1 = rot.'*rot1./nImages;
O2 = rot.'*rot2./nImages;

% We are registering one set of rotations (the estimated ones) to
% another set of rotations (the true ones). Thus, the transformation
% matrix between the two sets of rotations should be orthogonal. This
% matrix is either O1 if we recover the non-reflected solution, or O2,
% if we got the reflected one. In any case, one of them should be
% orthogonal.

err1 = norm(O1*O1.'-eye(3));
err2 = norm(O2*O2.'-eye(3));
% if (err1>TOL)
%     fprintf('Registering matrix is not orthogonal, err=%e  tol=%e\n',...
%         err1,TOL);
% end

errd1 = abs(abs(det(O1))-1); %allow O1 to be a reflection, i.e. in O(3)
errd2 = abs(abs(det(O2))-1); %allow O1 to be a reflection, i.e. in O(3)
% if (errd1>TOL)
%     fprintf('Determinant of registering matrix is not 1, err=%e  tol=%e\n',...
%         errd1,TOL);
% end

% In cany case, enforce the registering matrix O to be a rotation.
if err1<err2
    [U, ~, V]=svd(O1); % Use O1 as the registering matrix
    flag=1;
else
    [U, ~, V]=svd(O2); % Use O2 as the registering matrix
    flag=2;
end

O = V*U.';
rot_alligned = zeros(3,3,nImages);

% Plot estimation errors
diff=zeros(nImages,1);
for k=1:nImages
    Rref=q_to_rot(refq(:,k)).';
    
    R= g^sign_g_Ri(k).'*O*rotations(:,:,k);
    if flag ==2
        R = J*R*J;
    end
    rot_alligned(:,:,k) = R;
    
    diff(k) = norm(R - Rref,'fro');
end
mse= sum(diff.^2)/nImages;


return;

% if params.debug
%     figure,
%     hist(diff, 40);
%     set(gca,'FontSize', 20);
%     title('hist of err_{Fro} of reconstructed R_i')    
% end

end

function [view_dir_err_in_degrees,polar_ang_gt,stats_g] = check_viewdir_degrees_error(rots,refq)

g = [cosd(90) -sind(90) 0; 
	 sind(90)  cosd(90) 0; 
	 0 				 0  1]; % a rotation of 90 degrees about the z-axis
 
nImages = size(refq,2);
assert(size(rots,3) == nImages);

view_dir_err_in_degrees = zeros(1,nImages);
stats_g = zeros(1,nImages);
polar_ang_gt = zeros(1,nImages);
for i = 1:nImages
    
    Ri_gt  = q_to_rot(refq(:,i))';
    qi_gt = Ri_gt(:,3);
    
    polar_ang_gt(i) = acosd(qi_gt(3));
    
    qi = rots(:,3,i);
    qis = [qi,g*qi,g*g*qi,g*g*g*qi];
    
    
    ang_diffs_g = acosd(dot(repmat(qi_gt,1,4),qis));
    [YY,II] = min(ang_diffs_g);
    view_dir_err_in_degrees(i) = YY;
    stats_g(i) = II-1;
end
 
end

function [err_in_degrees,mean_err_per_rot] = check_degrees_error(rots,n_theta,refq)

FOUR = 4; % place-holder for a parameter that would be used for C3 as well
nImages = size(refq,2);
assert(size(rots,3) == nImages);
d_theta = 2*pi/n_theta;

g = [cosd(90) -sind(90) 0; 
	 sind(90)  cosd(90) 0; 
	 0 				 0  1]; % a rotation of 90 degrees about the z-axis
% precalculate g^s, s=0..params.FOUR-1
g_s = zeros(3,3,FOUR);
for s = 1:FOUR
    g_s(:,:,s) = g^(s-1);
end


err_in_degrees = zeros(1,n_theta*nImages);
mean_err_per_rot = zeros(1,nImages);
for k = 1:nImages
    rays_gt = zeros(n_theta,3);
    Rk_gt  = q_to_rot(refq(:,k))';
    for j = 1:n_theta
        rays_gt(j,:) = cos((j-1)*d_theta).*Rk_gt(:,1) ...
            +  ...
            sin((j-1)*d_theta).*Rk_gt(:,2);
    end
    
    c_err_in_degrees = cell(1,4);
    for s = 1:FOUR
        rays = zeros(n_theta,3);
        Rk   = g_s(:,:,s)*rots(:,:,k);
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
    Ri_gt = q_to_rot(refq(:,i)).';
    pol_ang(i) = acosd(Ri_gt(3,3));
end

% % err=acos(cos_angle);
% hist(err_in_degrees,50);
% h=gcf;
% set(gca,'FontSize',14);
% %title('Estimation error (in degrees)','FontSize',16);


end

function sign_g_Ri = g_sync(rotations, refq)

% The in-plain rotation stage retrieves all rotation matrices. However, every calculated rotation might 
% be a the rotated version of g^s, s=1,...,n from the ground truth. 
% Moreover, the rotation to be applied might be differernt for every image.
% We therefore synchronize all rotation so that only a possibly single global
% rotation should be applied to all rotation. Since the molecule is symmetric this doesn't 
%matter as far as the reconstruction is concerned, but it is important if 
% wish to compare to ground truth

nImages = size(rotations,3);
assert(nImages == size(refq,2));

g = [cosd(90) -sind(90) 0; 
	 sind(90)  cosd(90) 0; 
	 0 				 0  1]; % a rotation of 90 degrees about the z-axis

J = diag([1 1 -1]); % Reflection matrix

Ris_gt = zeros(3,3,nImages);
for i = 1:nImages
    Ris_gt(:,:,i) = q_to_rot(refq(:,i)).';
end

A_g = zeros(nImages,nImages);

for i=1:nImages 
    
    for j=i+1:nImages
        
        Ri  = rotations(:,:,i);
        Rj  = rotations(:,:,j);
        Rij = Ri.'*Rj;        
        
        Ri_gt = Ris_gt(:,:,i);
        Rj_gt = Ris_gt(:,:,j);
        diffs = zeros(1,4);
        for s=0:3
            Rij_gt = Ri_gt.'*(g^s)*Rj_gt;
            diffs(s+1) = min(norm(Rij-Rij_gt,'fro'),norm(Rij-J*Rij_gt*J,'fro'));
        end
        [~,ind] = min(diffs);
        A_g(i,j) = exp(-sqrt(-1)*2*pi/4*(ind-1));
    end
end
% A_g(k,l) is exp(-j(-theta_k+theta_j)) so we use transpose and not
% conjugate-transpose to obtain lower triangular entries
A_g = A_g + A_g';
% Diagonal elements correspond to exp(-i*0) so put 1. 
%This is important only for debug/verification purposes that spectrum is (K,0,0,0...,0)
A_g = A_g + eye(size(A_g)); 

[v,d] = eigs(A_g, 5, 'lm');
[evals,ind] = sort(diag(d), 'descend');
evect1 = v(:,ind(1));

log_message('Outer g syn first 5 eigenvalues=[%.2f %.2f %.2f %.2f %.2f]',...
    evals(1),evals(2),evals(3),evals(4),evals(5));

angles = exp(sqrt(-1)*2*pi/4*(0:3)).';
sign_g_Ri = zeros(nImages,1);
for ii  = 1:nImages
    zi  = evect1(ii);
    zi  = zi/abs(zi); % rescale so it lies on unit circle
    
    %Since a ccw and a cw closest are just as good, 
    % we take the absolute value of the angle
    angleDists = abs(angle(zi./angles));
    
    [~,I] = min(angleDists);
%     sign_g_Ri(ii) = 4-(I(1)-1);
    sign_g_Ri(ii) = I(1)-1;
end

end