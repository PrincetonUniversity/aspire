function [center, radgyr]=CenterOfMass(m)
% CENTEROFMASS Compute the center of mass of an image or volume
% 
% [center, radgyr]=CenterOfMass(m)
%   Find the center of mass and radius of gyration of a 2D (square) or 3D
%   (cube) density distribution.  The distribution is assumed to be
%   non-negative. 
%
% Based on code by Fred Sigworth, Yale University.
%
% Yoel Shkolnis,y November 2011.

dims=ndims(m);
m=m/sum(m(:));  % Scale so that it sums to unity.
n=size(m,1);    % All dimensions assumed to be the same

switch dims
    case 2
        if size(m,2)~=n
            error('2D input must be square');
        end
        
        [x,y]=ndgrid(-(n-1)/2:(n-1)/2);
        xc=x(:)'*m(:);
        yc=y(:)'*m(:);
        center=[xc yc];
        
        if nargout>1
            % Compute radius of gyration
            R2=(x-xc).^2+(y-yc).^2;
            radgyr=sqrt(R2(:)'*m(:));
        end;
        
    case 3
        
        if size(m,2)~=n || size(m,3)~=n
            error('3D input must be cube');
        end

        % First, compute center of mass
        [x, y, z]=ndgrid(-(n-1)/2:(n-1)/2);
        xc=x(:)'*m(:);
        yc=y(:)'*m(:);
        zc=z(:)'*m(:);
        center=[xc yc zc];
        
        if nargout>1
            % Compute radius of gyration
            R2=(x-xc).^2+(y-yc).^2+(z-zc).^2;
            radgyr=sqrt(R2(:)'*m(:));
        end;
end;