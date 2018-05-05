function shift_phases = calc_shift_phases(n_r,max_shift,shift_step)
%
%  Calculates the phases of all possible shifts
% 
% Input parameters:
%   n_r           Number of samples along each ray (in the radial direction
%   max_shift     (Optional) Maximal 1D shift (in pixels)  to search between
%                 common-lines. Default: 15.
%   shift_step    (Optional) Resolution of shift estimation in pixels. Default: 0.1
%
% Output parameters:
%   shift_phases  A matrix of size n_rxnshifts where shift_phases(:,shiftidx)
%                 is the fourier transform of corresponding shift 

if ~exist('shift_step','var')
    shift_step = 0.1;
end

if ~exist('max_shift','var')
    max_shift = 15;
end
              
% the maximum shift occurs at a diagonal direction
% max_shift = ceil(2*sqrt(2)*max_shift);
% there are max_shift in each direction
nshifts = 2*max_shift/shift_step+1;

rs = (0:n_r-1)'; % r_s should correspond to the real (underlying) frequency

shift_phases = zeros(n_r,nshifts);
for shiftidx = 1:nshifts
    shift = -max_shift + (shiftidx-1)*shift_step;
    % a shift in time domain corresponds to rotation in the frequency
    % domain
    shift_phases(:,shiftidx) = exp(-2*pi*sqrt(-1).*rs.*shift./(2*n_r-1));
end

end