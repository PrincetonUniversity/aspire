function PHI=register_orientations(PHI,refPHI)
%
% Register the estimated orienations to the reference orientations.
% The registration is required to eliminate the arbitrary rotation and
% reflection in the recovered orientations.
%
% PHI and refPHI are mx3 matrices.
%
% Yoel Shkolnisky, October 2007.

T = PHI\refPHI; % least squares to find unitary+reflection trasformation

% find best orthogonal approximation
[U,D,V] = svd(T);
Q = U * V';
PHI = PHI*Q; 
