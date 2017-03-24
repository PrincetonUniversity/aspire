function [euler_angles, r_valid_ks, r_good_ks, peak_width] =...
    cryo_sync3n_estimate_Rij(clstack, i, j, n_theta, VOTING_TICS_WIDTH, USE_MEDIAN_K)
% CRYO_SYNC3N_ESTIMATE_RIJ  Estimate relative rotation.
%
% cryo_sync3n_estimate_Rij(clstack, i, j, n_theta, VOTING_TICS_WIDTH, USE_MEDIAN_K)
%   Estimate relative rotations between image i and j. 
%   The function estimates relative rotation Ri'*Rj based on common lines
%   in clstack, where n_theta is the discretization of the common lines.
%   This function is called by cryo_sync3n_estimate_all_Rijs.
%   See cryo_sync3n_estimate_all_Rijs for a description of the parameters.
% 
% This function was renamed from  
%       R_pairs_parallel (clstack, i, j, algo)
%
% Yoel Shkolnisky, March 2016.
% Based on code by Ido Greenberg, 2015

% Initialization
N = size(clstack,1);
euler_angles = zeros(1, 3);

% First 2 Euler Angles induced directly by common lines.
euler_angles(3) = clstack(i,j);
euler_angles(1) = clstack(j,i);

% For the 3rd Euler Angle we need some other projection k.
% Look for valid k's, which have non singular geometry in relation to i,j:
phis = valid_k_for_ij(clstack, i, j, n_theta);
n_valid_ks = size(phis,1);

if n_valid_ks == 0   
    % No valid k's => define scores to vanish
    % GGG maybe try to prove projections are parallel and set some optimal
    % flat rotation between i,j - [* * 0 ; * * 0 ; 0 0 +-1]
    warning(['No valid k found for i,j=' int2str(i) ',' int2str(j) '. All triangles are too small. Projections might be parallel.']);
    euler_angles(2) = 1;
    r_valid_ks = 0;
    r_good_ks = 0;
    peak_width = 1;    
else    
    % Count valid k's
    r_valid_ks = n_valid_ks / (N-2);
    
    % Estimate angle between i,j planes
    [~, alpha, ~, width_needed_for_peak] = vote_k_for_ij(phis, VOTING_TICS_WIDTH);
    if USE_MEDIAN_K
        euler_angles(2) = median(alpha);
    else
        euler_angles(2) = mean(alpha);
    end
    
    % how wide peak was needed in order to find close angle estimations over k's,
    % and how many k's we found.
    r_good_ks = length(alpha) / (N-2);
    peak_width = width_needed_for_peak / VOTING_TICS_WIDTH;    
end
end
