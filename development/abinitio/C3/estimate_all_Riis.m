function Riis = estimate_all_Riis(self_cls,n_theta,refq)
%
% Finds the estimates for the self relative rotation matrices.
% Each estimate Rii is a rotation matrix that corresponds to either 
% R_i g R_i or to R_i g^3 R_i
% 
% 
% Input parameters:
%   self_cls   A 2-by-n table where n is the number of images. The i-th
%              column holds the indeces of the two self-common-lines in the
%              i-th image (with no particualr order)
%   n_theta    The angular resolution of each image. E.g., n_theta=360 means 
%              that each image has 360 lines 
%   refq       (Optional) A 4-by-n table. The i-th column represent the quaternion of
%              that corresponds to the rotation matrix of the i-th image
%
%
% Output parameters:
%   Riis       A 3-by-3-by-n array. The i-th slice is the 3-by-3 
%              rotation matrix Rii which is an estimate for either 
%              R_i g R_i or R_i g^3 R_i


if exist('refq','var') && ~isempty(refq)
    is_simulation = true;
else
    is_simulation = false;
end

log_message('Computing self relative orientations');

nImages = size(self_cls,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1: calculate the angle between Pi,Pgi (which is the same as the angle between Pi and Pg^3i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note : since we'll take the cosine of the angle we don't really care if
% we retrieve $-twiceAlpha$ or $360-twiceAlpha$ instead of $twiceAlpha$,
idxDiff      = self_cls(2,:)-self_cls(1,:);
% NOTE: DONT use method inds2clAngles. It assumes input to be positive (and
% even larger than zero), calculate therefore explitely the angle
twiceAlpha   = idxDiff*2*pi./n_theta;
cos_twiceAlpha = cos(twiceAlpha);

angles = acos(cos_twiceAlpha./(1-cos_twiceAlpha));
angles = real(angles); % TODO Gabi why this comes out complex???

% % we expect that cos_twiceAlpha lives in [-1,0], but sometimes due to
% % discretization issues we might get values larger than zero. Assert that
% % values are not too large and constrain them to be zero (these actually correspond to top-view images)
% assert(max(cos_twiceAlpha) <= eps,['max(cos_twiceAlpha) is ' num2str(max(cos_twiceAlpha))]);
% cos_twiceAlpha(cos_twiceAlpha>0) = 0;
% 
% angles = acos((1+cos_twiceAlpha)./(1-cos_twiceAlpha));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% debug code for step 1 (angles between planes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_simulation
    
    g = [cosd(120) -sind(120) 0; ...
     sind(120)  cosd(120) 0; ...
     0                 0  1]; % rotation matrix of 120 degress around z-axis

    assert(nImages == size(refq,2));
    
    angle_tol_err = 10/180*pi;
    
    angles_gt = zeros(1,nImages);
    for i=1:nImages
        Ri_gt = q_to_rot(refq(:,i))';
        Ri_3_gt = Ri_gt(:,3);
        angles_gt(i) = acos(Ri_3_gt.'*g*Ri_3_gt);
    end
    correct_idxs = find(abs(angles_gt-angles) < angle_tol_err);
    log_message('(P_i,gP_i) success rate=%.2f%%',numel(correct_idxs)/nImages*100);
    % find offending images (the tilt angle of their beaming direction)
    bad_idxs = find(abs(angles_gt-angles) >= angle_tol_err);
    bad_tilt_angles = zeros(1,numel(bad_idxs));
    k = 0;
    for bad_idx=bad_idxs
        k=k+1;
        Ri_bad = q_to_rot(refq(:,bad_idx))';
        bad_tilt_angles(k) = acosd(Ri_bad(3,3));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2: calculate the remaining euler angles moving from 
% between Pi,Pgi (which is the same as the angle between Pi and Pg^3i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aa = (self_cls(1,:)-1)*2*pi/n_theta;
bb = (self_cls(2,:)-1)*2*pi/n_theta + pi; % note the pi: C(gR_i,Ri) = C(R_i,g^3Ri)+pi.

% if ~params.real_data && strcmp(params.SCL,'GT') && params.confuse_scl_J
%     heads = find(round(rand(1,nImages)));
%     angles(heads) = -1*angles(heads);
% end

Riis = ang2orth(-bb, angles, aa);

for i=1:nImages
    [U,~,V] = svd(Riis(:,:,i)); %TODO: Gabi, is this needed? doesn't ang2orth return rotation matrices?
    Riis(:,:,i) = U*V.'; %Rii is the estimation of Ri^TgRi or Ri^Tg^3Ri up to J-conjugacy
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% debug code for step 2 (self relative orientations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if is_simulation
    
    J = diag([1 1 -1]); % Reflection matrix
    
    errs = zeros(1,nImages);
    errs_idx = zeros(1,nImages);
    for i=1:nImages
        Ri_gt = q_to_rot(refq(:,i))';
        Rii_g_gt  = Ri_gt.'*g*Ri_gt;
        Rii_g2_gt = Ri_gt.'*g^2*Ri_gt;
        
        Rii = Riis(:,:,i);
        [errs(i),errs_idx(i)] = min([norm(Rii-Rii_g_gt,'fro'),  norm(J*Rii*J-Rii_g_gt, 'fro'),...
            norm(Rii-Rii_g2_gt,'fro'), norm(J*Rii*J-Rii_g2_gt, 'fro')]);
    end
    cls_dist = histc(errs_idx,1:4)/numel(errs_idx);
    errs(bad_idxs) = [];
    log_message('MSE of Rii''s=%.2f',mean(errs.^2));
    log_message('cls_dist=[%.2f %.2f %.2f %.2f]',cls_dist);
    
%     [~,inds_err_Rii] = sort(errs,'descend');
    
%     params.inds_err_Rii = inds_err_Rii;
end

end