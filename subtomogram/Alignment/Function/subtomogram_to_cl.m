function [prop] = subtomogram_to_cl(subq, SNR, ECC)
% This function takes choices of parameters like SNR as inputs and output
% the estimated rotation between simulated subtomograms by applying
% selected methods.

% Inputs:
% subq   The quaternion representing the rotation between the subtomograms.
% SNR    The signal to noise ratio.
% alpha  The angular scale of the discretization of SO(3), the space of
% rotations in the sphere. Choices are from 6, 10, and 20 for now. Ideally
% 6 due to the radial resolution of 65.
% ECC    Choice of the score function for the common lines, where 0 stands
% for the squared euclidean distance as the MLE for Gaussian, and 1 stands
% for the cross-correlation between lines as their similarity.
% 
% Outputs:
% t      The time it takes for estimating the rotation between
% subtomograms, not including that for simulating the images.
% rot_est The estimated rotation matrix 
% rot_s   The true underlying rotation matrix
% err     The error between the estimated and the true rotation matrix,
% measured by the frobenius norm of their difference matrix.
% mi      The frobenius norm between the closest rotation matrix in the
% discretization of SO(3) to the true underlying rotation matrix, which is
% like the machine precision that cannot be overcome by any methods due to
% the resolution of such discretization.

% Yuan Liu, May 31, 2015

% STEP 0
% Generate simulated projections of size 65x65.
% For simplicity, the projections are centered.
n = 65;
ns = 9;
K_sub = 2;
K = K_sub*ns;

tiltq = zeros(4,K);   
grid = repmat((-4:4)*pi/24,1,K_sub);     % tilt angles in radians, their relative
                                        % angle is 7.5 deg = 1/24*pi
tiltq(1,:) = cos(grid/2);               
tiltq(2,:) = sin(grid/2);         % assume the rotation axis of tilt series 
                                % is the x-axis, the projection is in the z-axis
ref_q = zeros(4,K); % refq = b*a = q*r which is rotation b follows rotation a in absolute frame ### check
for i = 1:K
%     if i < 10
%         a = eye(3);
%     else
%         a = subq;
%     end
    a = subq(:,ceil(i/ns));
    b = tiltq(:,i);
    ref_q(1,i) = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4);
    ref_q(2,i) = a(1)*b(2) + a(2)*b(1) - a(3)*b(4) + a(4)*b(3);
    ref_q(3,i) = a(1)*b(3) + a(2)*b(4) + a(3)*b(1) - a(4)*b(2);
    ref_q(4,i) = a(1)*b(4) - a(2)*b(3) + a(3)*b(2) + a(4)*b(1);
end

[~,noisy_projs,~,~]=cryo_gen_projections_vsub(n,K,SNR,ref_q);
%masked_projs=mask_fuzzy(projs,33); % Applly circular mask
masked_noisy_projs = mask_fuzzy(noisy_projs,33);

n_theta=72;
n_r=65;
[npf,~]=cryo_pft(masked_noisy_projs,n_r,n_theta,'single');

max_shift=0;
shift_step=1;
[ref_clstack,~]=clmatrix_cheat_q(ref_q,n_theta);
tol = 360/n_theta;
if ECC == 1
    [clstack, ~, ~, ~] = commonlines_gaussian_vsub(npf,max_shift,shift_step);
    prop=comparecl( clstack, ref_clstack, n_theta, tol );
else
    [clstack_e, ~, ~, ~] = commonlines_gaussian_vsubeucl(npf,max_shift,shift_step);
    prop=comparecl( clstack_e, ref_clstack, n_theta, tol );
end


end

