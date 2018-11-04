function idxs = clAngles2Ind(clAngles,n_theta)
%
% Transforms pairs [x,y] of spatial corrdinates into angle indexes
% 1,2,...,n_theta representing angles [0,2pi]
% 
% Input parameters:
%   clAngles   A 2xn table where each column [x,y] represents the patial
%              coordinate of some point
%   n_theta    Angular resolution. Number of Fourier rays computed for each
%              projection.
%
% Output parameters:
%   idxs       An array of length n. The i-th entry is a number inn [1,n_theta] 
%              representing the angle index of clAngles(:,i) 

assert(size(clAngles,1) == 2);

thetas = atan2(clAngles(2,:),clAngles(1,:));

thetas = mod(thetas,2*pi); % Shift from [-pi,pi] to [0,2*pi).

idxs = thetas./(2*pi)*n_theta; % linear scale from [0,2*pi) to [0,n_theta).
idxs = mod(round(idxs),n_theta)+1;

end