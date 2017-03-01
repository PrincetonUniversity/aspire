function [ R , r_valid_ks , r_good_ks , peak_width , euler_angles ] =...
    cryo_sync3n_estimate_all_Rijs(clstack, n_theta, VOTING_TICS_WIDTH,USE_MEDIAN_K,progress_update_freq)
% CRYO_SYNC3N_ESTIMATE_ALL_RIJS  Estimate all relative rotations.
%
% cryo_sync3n_estimate_all_Rijs(clstack, n_theta, VOTING_TICS_WIDTH,USE_MEDIAN_K,progress_update_freq)
% Estimate all relative rotations Ri'*Rj based on common lines matrix
% clstack. The angular resolution of the common lines is n_theta.
% 
% Input:
%   clstack Common lines matrix (clstack(i,j) = the index of the common
%           line between i&j in the image of i)
%   n_theta Discretization of common lines
%   VOTING_TICS_WIDTH   Resolution of angles (in degrees) used in k's
%           voting algorithm for the angle between projections. Default is
%           1 degree.
%   USE_MEDIAN_K    When estimating the angle between projections (i,j)
%           using all third images k, whether to use median or mean
%           estimations over {k}. Default is 0 (use mean).
%   progress_update_freq    Frequency of progress update prints (in loop
%           iterations). Default is 25000.
% 
% Output:
%   R   3X3XN_pairs array of estimated relative rotations. 
%   r_valid_ks  For every pair (i,j), the number of third images k which have
%               non-singular geometry in relation to (i,j).
%   peak_width 	For every pair (i,j), how wide peak was needed to include
%               entries in the (i,j)-angle-estimations-histogram. 
%               XXX I don't understand this description. XXX
%   r_good_ks   For every pair (i,j), the number of thirf images k in the
%               peak of the histogram.
%   euler_angles    For every pair (i,j), the 3 Euler angles representing
%   the relative rotation between i&j.
% 
% This function was renamed from 
%       R_pairs(algo, clstack, progress_update_freq)
%
% Yoel Shkolnisky, March 2016
% Baed on code by Ido Greenberg, 2015


if ~exist('VOTING_TICS_WIDTH','var') || isempty(VOTING_TICS_WIDTH)
    % Voting histogram resolution (in degrees).
    % The sensitivity to this variable was not examined. Historically, the values 1 and 3 were in use.
    % The typical errors in the common lines are known to be in the range of 5-20 degrees, so 3 degrees seems like a fair default resolution for the histogram.
    VOTING_TICS_WIDTH=3;
end
log_message('VOTING_TICS_WIDTH=%d',VOTING_TICS_WIDTH);


if ~exist('USE_MEDIAN_K','var') || isempty(USE_MEDIAN_K)
    USE_MEDIAN_K=0;
end
log_message('USE_MEDIAN_K=%d',USE_MEDIAN_K);

if ~exist('progress_update_freq','var') || isempty(progress_update_freq)
    progress_update_freq=25000;
end
log_message('progress_update_freq=%d',progress_update_freq);

% Initialization
progress_update_next = progress_update_freq;
N = size(clstack,1);
N_pairs = N*(N-1)/2;

euler_angles = zeros(N_pairs, 3);
r_valid_ks = zeros(N_pairs, 1);
r_good_ks = zeros(N_pairs, 1);
peak_width = zeros(N_pairs, 1);

% Compute relative rotations
idx = 0;
for i = 1:(N-1)
    
    % INEFFICIENT IMPLEMENTATION: The parfor loop probably can be wrapped in more elegant and efficient way (as in similar parallel loops in previous variants of the algorithm).
    % However, it does not expected to affect the running time significantly.
    
    euler_angles_tmp = cell(N-i,1);
    r_valid_ks_tmp = cell(N-i,1);
    r_good_ks_tmp = cell(N-i,1);
    peak_width_tmp = cell(N-i,1);

    % Calculate
    parfor j = (i+1):N
        [euler_angles_tmp{j-i}, r_valid_ks_tmp{j-i}, r_good_ks_tmp{j-i}, peak_width_tmp{j-i}] =...
            cryo_sync3n_estimate_Rij(clstack, i, j, n_theta,VOTING_TICS_WIDTH,USE_MEDIAN_K)
    end

    % Sync calculations
    for j = (i+1):N
        idx = idx + 1;
        euler_angles(idx, :) = euler_angles_tmp{j-i};
        r_valid_ks(idx) = r_valid_ks_tmp{j-i};
        r_good_ks(idx) = r_good_ks_tmp{j-i};
        peak_width(idx) = peak_width_tmp{j-i};
    end

    % Print updates
    if progress_update_freq > 0 && idx >= progress_update_next
        log_message('Computed %d Rij''s out of %d', progress_update_next, N_pairs);
        while progress_update_next <= idx
            progress_update_next = progress_update_next + progress_update_freq;
        end
    end

end

% Convert common lines to Euler angles
euler_angles(:, 3) = (euler_angles(:,3)-1) * 2*pi/n_theta - pi;
euler_angles(:, 1) = pi - (euler_angles(:,1)-1) * 2*pi/n_theta;

% Euler Angles --> Rotations
R = ang2orth( euler_angles(:,1) , euler_angles(:,2) , euler_angles(:,3) );

end
