% Isotropic two dimnetional window 
%
% Syntax: [ w ] = iwindow( p, 'hann' )
%
% Inputs:
%    p - size of window
%    kind - window kind {'hamming', 'hann'}
%    win_params - optional window parameters
%
% Outputs:
%    w - window samples 
%
%
% Other m-files required: cart2rad.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Cohen, Mor
%________________________________________

% email: morcohe5@mail.tau.ac.il
% May 2011; Last revision: 04-May-2011
%
% See also: 

function [ w ] = iwindow( p, kind )
I = cart2rad(p);
M = p/2;
switch kind
    case 'hamming'
        w = 0.54-0.46.*cos(2*pi*(I+M)./p);
    case 'hann'
        w = 0.5*(1-cos(2*pi*(I+M)./p));
    case 'boxcar'
        w = ones(p);
    case 'bartlett'
        w=1-abs(I)/floor((p-1)/2);
        w=max(w,0);
    otherwise
        error('No such window %s',kind);
end
w(I>(M+1)) = 0;
end

