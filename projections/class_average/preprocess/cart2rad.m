% Compute centralized image coordinates in polar radius size values
%
% Syntax: [ I ] = cart2rad(N)
%
% Inputs:
%    N - the size of two dimentional lattice NxN
%
%
% Outputs:
%   I - centralized polar samples
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:  Mor Cohen
%________________________________________

% email: morcohe5@mail.tau.ac.il
% May 2011; Last revision: 04-May-2011
%
% See also: 
function [ I ] = cart2rad(N)
    N = floor(N);    
    p = (N-1)/2;
    [X,Y] = meshgrid(-p:p,-p:p);
    
    I = sqrt(X.^2+Y.^2);
end

