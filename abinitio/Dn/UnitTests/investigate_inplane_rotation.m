function [Ris_tilde, theta_est, original_ests] = investigate_inplane_rotation(npf,vis,...
    n_symm, inplane_rot_res)
%
% General description
% 
% Input parameters:
%   npf              A 3D array where each image npf(:,:,i) corresponds to the Fourier
%                    transform of projection i.
%   vis              description
%   inplane_rot_res  (Optional) The resolution in degrees of in-plane rotation angle.
%                    E.g., resolution=1 means we look for angles 0,1,2,3 ...
%                    E.g.,  resolution=0.25 means we look for angles 0,0.25,0.5,...
%                    Defualt=1
%   max_shift 
%   shift_step          
% Output parameters:
%   rots       A 3D array of size 3x3xnImages where the i-th slice
%              rots(:,:,3) is equal to Ri (the i-th rotation matrix)

if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('max_shift','var')
    max_shift = 0;
end

if ~exist('inplane_rot_res','var')
    inplane_rot_res = 1;
end


log_message('estimating in-plane rotation angles')

[n_r,n_theta,nImages] = size(npf);

% create all in-plane angular differences using 'inplane_rot_res' as the
% discretization

max_angle = floor(360/n_symm)*n_symm;

theta_ij = (0:inplane_rot_res:(max_angle-inplane_rot_res))*pi/180;
n_theta_ij = numel(theta_ij);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: construct all in-plane rotation matrices R(theta_ij)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gx = diag([1, -1, -1]);
R_theta_ij = quat2rotm([cos(theta_ij / 2)', zeros(n_theta_ij, 2), ...
    sin(theta_ij / 2)']);

counter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: construct all tilde_Ri. These are rotation matrices 
% whose third row is equal to that of Ri)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ris_tilde = zeros(3,3,nImages);
for i = 1:nImages
    v_i = vis(:,i).';
    Ris_tilde(:,:,i) = complete_3rdRow_to_rot(v_i);
end
original_ests = Ris_tilde;

% ignoring dc term
npf_tmp = npf;
npf_tmp(1,:,:) = 0;

npf_normalized = bsxfun(@rdivide,npf_tmp,...
                                 sqrt(sum((abs(npf_tmp)).^2)));

clear npf_tmp;

%precompile the shift phases
shift_phases = calc_shift_phases(n_r,max_shift,shift_step);
[~,nshifts] = size(shift_phases);
npf_i_shifted = zeros(n_r,n_theta,nshifts);
estimated_angles = cell(1, nImages);
estimations_mat = NaN(nImages, nImages, 2);
tol = 1e-12;
angle_mod = 2 * pi / n_symm;
angle_mod_deg = rad2deg(angle_mod);
max_angle = angle_mod - deg2rad(inplane_rot_res);
%max_angle_deg = rad2deg(max_angle);
max_index = find(are_almost_equal(theta_ij, max_angle, tol));
angles_degrees = rad2deg(theta_ij);
all_scores_table = NaN(nImages, nImages, max_index, max_index);
theta_diff = zeros(nImages);
theta_sum = zeros(nImages);

for i = 1:nImages
    
    npf_i = npf(:,:,i);
    % get all possible shifted copies of the image
    for s=1:nshifts
        npf_i_shifted(:,:,s) = bsxfun(@times,npf_i,shift_phases(:,s));
    end
    
    npf_i_shifted(1,:) = 0; % ignoring dc term
    % normalize each ray to be of norm 1
    npf_i_shifted = bsxfun(@rdivide,npf_i_shifted,...
        sqrt(sum((abs(npf_i_shifted)).^2)));
    
    Npf_i_shifted = gpuArray(single(npf_i_shifted));
    
    % for performance purposes we precalculate this once per each image i.
    % If it were for readibility it would be better to put this in the inner
    % loop.
    Ri_tilde_t = Ris_tilde(:, :, i).';
    % calculate R_i_tilde*R_theta_ij for all possible R_theta_ijs
    tmp = multiprod(Ri_tilde_t, R_theta_ij);
    tmp_gx = multiprod(Ri_tilde_t * gx, R_theta_ij);
    
    for j = 1:nImages
        
%         t1 = clock;
        counter = counter+1;
        
        npf_j = npf_normalized(:,:,j);        
        
        Rj_tilde = Ris_tilde(:, :, j);
        
        % to each of the previous computed matrices we multiply it by Rj_tilde
        Us  = multiprod(tmp,Rj_tilde);
        Us_gx = multiprod(tmp_gx, Rj_tilde);
        
        Us = reshape(Us, 9, n_theta_ij);
        Us_gx = reshape(Us_gx, 9, n_theta_ij);
        % extract a common-line index for each possible theta_ij
        c_1 = [-Us(8,:) ;  Us(7,:)]; %[-U(2,3)  U(1,3)];
        c_2 = [ Us(6,:) ; -Us(3,:)]; %[ U(3,2) -U(3,1)];
        c_1_gx = [-Us_gx(8,:) ;  Us_gx(7,:)]; %[-U(2,3)  U(1,3)];
        c_2_gx = [ Us_gx(6,:) ; -Us_gx(3,:)]; %[ U(3,2) -U(3,1)];
        
        cijs = clAngles2Ind(c_1,n_theta);
        cjis = clAngles2Ind(c_2,n_theta);
        cijs_gx = clAngles2Ind(c_1_gx,n_theta);
        cjis_gx = clAngles2Ind(c_2_gx,n_theta);
        Npf_j = gpuArray(single(npf_j));
        
        % cross-correaltion - so we want to conjugate
        co = bsxfun(@times,Npf_i_shifted(:,cijs,:),conj(Npf_j(:,cjis))); 
        Corrs = sum(co);
        corrs = gather(Corrs);        
        corrs = reshape(corrs, n_theta_ij/n_symm,n_symm,nshifts);
        
        co_gx = bsxfun(@times,Npf_i_shifted(:,cijs_gx,:),conj(Npf_j(:,cjis_gx))); 
        Corrs_gx = sum(co_gx);
        corrs_gx = gather(Corrs_gx);        
        corrs_gx = reshape(corrs_gx, n_theta_ij/n_symm,n_symm,nshifts);
        
        
        if nshifts > 1
            % the 1d shifts depend not only on the 2d shifts of each and every
            % image, but also on their common-lines. That is, for a given image, the 1d
            % shift might be different for different lines !
            % In particular, this means that each one of the 4 lines we are considering in a given image
            % may correspond to a DIFFERENT 1d shift. Therefore each line has the dof, so to speak, of choosing its own shift.
            % All in all, this means that we max out the shift BEFORE taking
            % the mean score over all n lines.
            corrs = max(real(corrs),[],3);
            corrs_gx = max(real(corrs_gx), [], 3);
        end
        % now take the mean score over all groups of n-pairs-of-lines , and find the group that attains the maximum
        corrs = mean(real(corrs), 2);
        corrs_gx = mean(real(corrs_gx), 2);
        
        scores_table = correlation_scores(angles_degrees, corrs, ...
            corrs_gx, max_index, angle_mod_deg, n_symm, tol);
        [max_score, I] = max(scores_table(:));
        [theta_i_ind, theta_j_ind] = ind2sub([max_index, max_index], I);
        theta_i = angles_degrees(theta_i_ind);
        theta_j = angles_degrees(theta_j_ind);
        all_scores_table(i, j, :, :) = scores_table;
        theta_diff(i, j) = -theta_i + theta_j;
        theta_sum(i, j) = theta_i + theta_j;
%         theta_diff(j, i) = -theta_diff(i, j);
%         theta_diff(i, i) = 0;
    end
end

%theta_sum = (theta_sum + theta_sum.') / 2;
%theta_diff = (theta_diff - theta_diff.') / 2;
cosine_mat = (cosd(n_symm * theta_diff) + cosd(n_symm * theta_sum)) / 2;
cosine_mat = (cosine_mat + cosine_mat.') / 2;
sine_mat = (sind(n_symm * theta_diff) + sind(n_symm * theta_sum)) / 2;
[v, d] = eigs(cosine_mat, 7);
[d, idx] = sort(diag(d), 'descend');
%d = diag(d);
disp(d);
leading_eigenvalue = d(1);
cosine_vec = sqrt(leading_eigenvalue) * v(:, idx(1));
cosine_vec = min(max(cosine_vec, -1), 1);
[U, S, V] = svds(sine_mat, 7);
%sine_mat = U(:, 1) * S(1, 1) * V(:, 1)';
sine_vec = (cosine_vec.' * sine_mat).' / leading_eigenvalue;
sine_vec = min(max(sine_vec, -1), 1);
sine_vec2 = estimate_sine_vec(sine_mat, cosine_vec, nImages, leading_eigenvalue);
fprintf('Diff sine:\n');
disp(norm(sine_vec - sine_vec2));
%sine_vec = V(:, 1) * S(1, 1) / sqrt(leading_eigenvalue);
disp(diag(S));
% sine_vec_alternative = sqrt(1 - cosine_vec.^2);
% replaced_indices = abs(sine_vec - sine_vec_alternative) < ...
%     abs(sine_vec + sine_vec_alternative);
% sine_vec(replaced_indices) = sine_vec_alternative(replaced_indices);
theta_est = rad2deg(atan2(sine_vec, cosine_vec)).' / n_symm;

log_message('forming rotation matrices');
in_plane_rots = quat2rotm(quaternion([zeros(nImages, 2), theta_est.'], ...
   'eulerd', 'XYZ', 'frame'));
Ris_tilde = multiprod(in_plane_rots, Ris_tilde, [1 2], [1 2]);
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

function scores = correlation_scores(angles_list, corrs, corrs_gx, ...
    max_index, angle_mod, n_symm, tol)
    %% Computation of the scores table for all angles' pairs.

    scores = NaN(max_index, max_index);
    for i = 1:max_index
        for j = 1:max_index
            % Picking two angles.
            theta_i = angles_list(i);
            theta_j = angles_list(j);
            
            % Computing the score for these specific angles.
            total_score = common_lines_score(angles_list, theta_i,...
                theta_j, corrs, corrs_gx, angle_mod, tol);
            scores(i, j) = total_score;
            
            % Validating the main identity of the score function.
            assert_adjoint_angles(theta_i, theta_j, corrs,...
                corrs_gx, angles_list, n_symm, total_score, ...
                angle_mod, tol, 0);
        end
    end
end

function score = common_lines_score(angles_vec, theta_i, theta_j, ...
    corrs, corrs_gx, angle_mod, tol)
    %% Evaluation of the common-lines score function in two-variables.
    
    % Finding the indices of the sum and the difference of the angles.
    sum_index = find_angle_mod(angles_vec, theta_i + theta_j, ...
        angle_mod, tol);
    diff_index = find_angle_mod(angles_vec, -theta_i + theta_j, ...
        angle_mod, tol);
    
    % Computing the score.
    score = max([corrs(sum_index) + corrs_gx(diff_index),...
                corrs(diff_index) + corrs_gx(sum_index)]);
end

function score = common_lines_single_param_score(angles_list, theta_i, ...
    corrs, corrs_gx, angle_mod)

    sum_index = find_angle_mod(angles_list, 2 * theta_i, angle_mod);        
    score = max([corrs_gx(sum_index), corrs(sum_index)]);
end

function sine_vec = estimate_sine_vec(sine_mat, cosine_vec, nImages, lead_eigen)
x0 = zeros(nImages, 1);
A = [];
b = [];
Aeq = [];
beq = [];
lb = -1;
ub = 1;

function [f, grad] = sine_mat_fit(x)
    f = norm(cosine_vec * x' - sine_mat, 'fro') / 2;
    grad = lead_eigen * x - sine_mat' * cosine_vec;
end
hess_func = @(x, lambda) lead_eigen * eye(nImages);

function [c, cesq, gradc, gradcesq] = unit_circle_contr(x)
    cesq = [];
    gradcesq = [];
    c = x.^2 + cosine_vec.^2 - 1;
    gradc = diag(2 * x);
end

options = optimoptions('fmincon','Algorithm','interior-point',...
    "SpecifyConstraintGradient",true,"SpecifyObjectiveGradient",true...
    , 'HessianFcn',hess_func);

fun = @sine_mat_fit;
disp(fun(x0));
[sine_vec,fval,exitflag,output, grad, h] = fmincon(@sine_mat_fit,x0,A,b,...
    Aeq,beq,lb,ub,@unit_circle_contr,options);
end

function scores = self_common_line_scores(angles_list, corrs, corrs_gx, ...
    max_index, angle_mod)

    scores = NaN(1, max_index);
    for i = 1:max_index
        theta_i = angles_list(i);        
        scores(i) = common_lines_single_param_score(angles_list, theta_i, ...
            corrs, corrs_gx, angle_mod);
    end
end

function assert_score_function_properties(scores_table, max_score, ...
    theta_i, n_symm, is_single_param_score)

    % Verify that there are exactly four optimal pairs of angles
    % for the common-lines score function.
    tol = 1e-14;
    optimal_points_count = numel(find(are_almost_equal(...
        scores_table, max_score, tol)));
    if is_single_param_score
        assert(optimal_points_count == 2 || optimal_points_count > 70);
    elseif are_almost_equal(theta_i, 0, tol) || ...
            are_almost_equal(theta_i, 180 / n_symm, tol)
        assert(optimal_points_count == 2);
    else
        assert(optimal_points_count == 4);
    end    
end

function assert_adjoint_angles(theta_i, theta_j, corrs, corrs_gx, ...
    angles_degrees, n_symm, score, angle_mod, tol, is_single_param_score)
    %% This function verifies the equality of score between (theta_i, theta_j)
    %% and (pi/n + theta_i, pi/n + theta_j) for given angles
    adjoint_theta_i = mod(180 / n_symm + theta_i, angle_mod);
    adjoint_theta_j = mod(180 / n_symm + theta_j, angle_mod);
    adjoint2_theta_i = mod(360 / n_symm - theta_i, angle_mod);
    adjoint3_theta_i = mod(180 / n_symm - theta_i, angle_mod);
    
    if ~is_single_param_score
        total_score = common_lines_score(angles_degrees, adjoint_theta_i,...
            adjoint_theta_j, corrs, corrs_gx, angle_mod, tol);
        total_score2 = common_lines_score(angles_degrees, adjoint2_theta_i,...
            theta_j, corrs, corrs_gx, angle_mod, tol);
        total_score3 = common_lines_score(angles_degrees, adjoint3_theta_i,...
            adjoint_theta_j, corrs, corrs_gx, angle_mod, tol);

        assert(abs(score - total_score2) < 1e-13);
        assert(abs(score - total_score3) < 1e-13);
    else
        total_score = common_lines_single_param_score(angles_degrees, ...
            theta_i, corrs, corrs_gx, angle_mod); %max([corrs_gx(sum_index), corrs(sum_index)]);
    end
    assert(abs(score - total_score) < 1e-13);
end

function assert_adjoint_angles2(theta_i, theta_j, corrs, corrs_gx, ...
    angles_degrees, n_symm, score, angle_mod, tol, is_single_param_score)
    %% This function verifies the equality of score between (theta_i, theta_j)
    %% and (pi/n + theta_i, pi/n + theta_j) for given angles
    adjoint_theta_i = mod(180 / n_symm + theta_i, angle_mod);
    adjoint_theta_j = mod(180 / n_symm + theta_j, angle_mod);
    adjoint2_theta_j = mod(360 / n_symm - theta_j, angle_mod);
    adjoint3_theta_j = mod(180 / n_symm - theta_j, angle_mod);
    
    if ~is_single_param_score
        total_score = common_lines_score(angles_degrees, adjoint_theta_i,...
            adjoint_theta_j, corrs, corrs_gx, angle_mod, tol);
        total_score2 = common_lines_score(angles_degrees, theta_i,...
            adjoint2_theta_j, corrs, corrs_gx, angle_mod, tol);
        total_score3 = common_lines_score(angles_degrees, adjoint_theta_i,...
            adjoint3_theta_j, corrs, corrs_gx, angle_mod, tol);

        assert(abs(score - total_score2) < 1e-13);
        assert(abs(score - total_score3) < 1e-13);
    else
        total_score = common_lines_single_param_score(angles_degrees, ...
            theta_i, corrs, corrs_gx, angle_mod); %max([corrs_gx(sum_index), corrs(sum_index)]);
    end
    assert(abs(score - total_score) < 1e-13);
end

function result = are_almost_equal(a, b, err)
    result = abs(a - b) <= err;
end
    