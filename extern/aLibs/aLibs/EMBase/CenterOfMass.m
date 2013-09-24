function [center, radgyr]=CenterOfMass(m)
% function [center, radgyr]=CenterOfMass(m)
% Find the center of mass and radius of gyration of a 2D (square) or 3D
% (cube) density distribution.  The distribution is assumed to be non-negative.

dims=ndims(m);
m=m/sum(m(:));  % scale so that it sums to unity.
n=size(m,1);    % all dimensions assumed to be the same

switch dims
    case 2
        [x,y]=ndgrid(floor(-(n-1)/2):floor((n-1)/2));
        xc=x(:)'*m(:);
        yc=y(:)'*m(:);
        center=[xc yc];
        
        if nargout>1
            % Compute radius of gyration
            R2=(x-xc).^2+(y-yc).^2;
            radgyr=sqrt(R2(:)'*m(:));
        end;
        
    case 3
        % First, compute center of mass
        [x, y, z]=ndgrid(floor(-(n-1)/2):floor((n-1)/2));
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