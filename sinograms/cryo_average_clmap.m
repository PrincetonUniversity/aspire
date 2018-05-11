function avg_map = cryo_average_clmap(map, d, W0)
%
% CRYO_AVERAGE_CLMAP    Apply filtering to improve common line detection
%
% cryo_average_clmap(map, d)    Filter the correlation matrix of a pair of
%   (Fourier transformed projections) by averaging each correaltion with
%   its neighbors in a dxd neighborhood. This is based on the idea that not
%   only the pair common lines should have high correlations, but also
%   their neighboring lines. Defaults weights W0 are used.
%
% cryo_average_clmap(map, d, W0)    Use the given W0 as the weights for the
%   weighted average.
%
% map is an Lx2L matrix, and the functions averagea every 2d+1 X 2d+1
% square of it, using weighted average defined by W0.
%
% Ido Greenberg 2016.

% input validation
if ~(size(map,2)==2*size(map,1)); error('invalid map size, must be L on 2L'); end
if ~(d>0); error('invalid square size, must be d>0'); end

% configure weights
if ~exist('W0','var'); W0 = 1.25.^(21:-1:1); end
W0 = W0(1:(1+2*d));
W = W0 ( 1 + kron(abs(-d:d),ones(2*d+1,1)) + kron(abs(-d:d)',ones(1,2*d+1)) );

% extend map
L = size(map,1);
map = [map ; map(:,[L+1:2*L,1:L])];

% average map
avg_map = ifft2(...
    fft2(map) .* conj(fft2(W,size(map,1),size(map,2))) ...
    );

% Reshape map.
% Currently, every pixel is the average of a square, where the pixel is
% the upper-left corner of the square. We want to shift the pixels by d
% to the center of the corresponding squares.
% this turns into:
%   [ avg_map(2L+1-d:2L,2L+1-d:2L) , avg_map(2L+1-d:2L,1:2L-d) ;
%   avg_map(1:2L-d,2L+1-d:2L) , avg_map(1:2L-d,1:2L-d) ]
%   or equivalently:
% avg_map( [ 2L+1-d:2L , 1:2L-d ] , [ 2L+1-d:2L , 1:2L-d ] ).
% We only need Lx2L of of the current LxL map. This turns into
% avg_map(1:L,1:2*L). together, we have:
avg_map = avg_map([2*L+1-d:2*L,1:L-d],[2*L+1-d:2*L,1:2*L-d]);

end
