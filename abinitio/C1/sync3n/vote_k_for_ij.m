
function [ good_ks, alpha , peakh , w_theta_needed ] = vote_k_for_ij( phis , w_theta )
%VOTE_k_FOR_ij Compute the histogram of the angles between
% projections i and j, usings different k's.
% Pick the angles closest to the peak.
%
% input:
% phis - array of cosines of angles induced by valid k's, and the valid k's indices
% w_theta - width of each angle tic [degrees]
% 
% output:
% good_ks - number of k's in the peak of the histigram
% alpha - estimated angles near the peak
% peakh - maximum of the angle estimations histogram
% w_theta_needed - how wide peak was needed in order to find close angle estimations

% Parameters used to compute the smoothed angle histogram. Same as in 
% histfilter_v6 (for compatibility).
angle_tics = 0:w_theta:180;
h = zeros(size(angle_tics));

angles = acos(phis(:,1))*180/pi;
% Angles that are up to 10 degrees apart are considered
% similar. This sigma ensures that the width of the density
% estimation kernel is roughly 10 degrees. For 15 degress, the
% value of the kernel is negligible.

angles_distances = zeros(size(phis,1),numel(angle_tics));
for i = 1:numel(angle_tics)
    % how close are the angles to the current tic?
    angles_distances(:,i) = angle_tics(i)-angles;
    %h(i) = gather(sum(exp(-(angle_tics(i)-angles).^2/(2*sigma.^2))));
end

%sigma=2.64;
sigma = 3.0; % For compatibility with histfilter_v6
h = sum(exp(-angles_distances.^2/(2*sigma.^2)));
%h = gather(sum(exp(-gpuArray(angles_distances).^2/(2*sigma.^2))));
% GPU was found to be not cost-effective in this case - wastes more time than it saves.
% Alternative implementation can use bsxfun to compute h=exp(angles*tics-angles.^2-tics.^2) (which is algebraically equivalent). We did not check which implementation is more efficient.
% This is quite important since this specific line is a significant bottle-neck in the running time.

% We assume that at the location of the peak we get the true angle
% between images k1 and k2. Find all third images k3, that induce an
% angle between k1 and k2 that is at most 10 off the true
% angle. Even for debugging, don't put a value that is smaller than two
% tics, since the peak might move a little bit due to wrong k3 images
% that accidentally fall near the peak.
[peakh, peakidx] = max(h);

% look for the estimations in the peak of the histogram
w_theta_needed = 0;
idx = [];
while numel(idx) == 0
    w_theta_needed = w_theta_needed + w_theta; % widen peak as needed
    idx = find( abs(angles-angle_tics(peakidx)) < w_theta_needed );
end

% k's in the peak
good_ks = phis(idx,2);

% Good angles [radians]
alpha = acos(phis(idx,1));

end
