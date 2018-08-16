% PFT_FREQS Polar Fourier transform frequencies
%
% Usage
%    freqs = pft_freqs(r, n_theta);
%
% Input
%    r: An array of radial frequencies.
%    n_theta: The number of samples on each concentric circle. Must be even.
%
% Output
%    freqs: An array of size n_r*n_theta/2, where n_r is the length of r,
%       containing the coordinates of the frequency samples between the angles
%       of 0 and pi.

% Written by Zhizhen Zhao - 02/2015
% Reformatted and documented by Joakim Anden - 2018-Apr-13
% Optimized by Joakim Anden - 2018-Jul-24

function freqs = pft_freqs(r, n_theta)
	if mod(n_theta, 2) ~= 0
		error('`n_theta` must be even');
	end

    n_r = length(r);
    dtheta = 2*pi/n_theta;

    % Sampling points in the Fourier domain
    freqs = zeros(n_r*n_theta/2, 2);
    for j = 1:n_theta/2
        freqs((j-1)*n_r+1:j*n_r,:) = ...
            [r*sin((j-1)*dtheta), r*cos((j-1)*dtheta)];
    end
end
