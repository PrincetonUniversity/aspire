function vol=cryo_gaussian_phantom_3d(def,n,rmax)
%
% 3D Guassian phantom.
% The phantom is a linear combination of 3D Guassians.
%
% Input parameters:
%   def     Filename of the efinitions file of the phantom.
%   n       Size of the output phantom (containing n^3 volxels)
%   rmax    Physical dimensions of the phantom. The phantom has n samples
%           in each dimension between -rmax and rmax. 
%           If this functions is combined with cryo_project_gaussian, then
%           rmax must match rmax used in cryo_project_gaussian for
%           generating proejctions.  
% 
% Output parameters:
%   vol     The 3D phantom. An nxnxn volume, where each voxel is equal to
%           the sum of the values of all Gaussians at that voxel.
%  
% Examples:
%
%   vol = cryo_gaussian_phantom_3d('C1_params',128,1);
%   figure; imshow(squeeze(vol(64,:,:)));
%
% Yoel Shkolnisky, June 2007.
%
% Y.S May 2013 Change line 39 from 
%                       gparams=gaussian_params3;
%              to       gparams=C1_params;
%            
%               Revise input format.
%
% Y.S. August 2013 Cosmetic changes


% param can be either a numeric scalar of a 3 by m array
if (nargin~=3) || ~ischar(def) || ~isscalar(n)
    error('Input must be (def,n,rmax)');
end


gparams=feval(def); %load the approriate phantom.

vol = zeros([n n n]);

if mod(n,2)==1
    rng=(-(n-1)/2:(n-1)/2) / ((n-1)/2)*rmax; 
else
    rng=(-n/2+1/2:n/2-1/2)/(n/2)*rmax;
end

[x,y,z] = ndgrid(rng,rng,rng); 

coord = [flatten(x); flatten(y); flatten(z)];
vol = flatten(vol);

for k = 1:size(gparams,1)    
   rho = gparams(k,1);          % Amplitude change for this ellipsoid
   asq = gparams(k,2);          % std x
   bsq = gparams(k,3);          % std y
   csq = gparams(k,4);          % std z
   x0 = gparams(k,5);           % x offset
   y0 = gparams(k,6);           % y offset
   z0 = gparams(k,7);           % z offset
   phi = gparams(k,8)*pi/180;   % first Euler angle in radians
   theta = gparams(k,9)*pi/180; % second Euler angle in radians
   psi = gparams(k,10)*pi/180;  % third Euler angle in radians
   
   cphi = cos(phi);
   sphi = sin(phi);
   ctheta = cos(theta);
   stheta = sin(theta);
   cpsi = cos(psi);
   spsi = sin(psi);
   
   % Euler rotation matrix
   A = [cpsi*cphi-ctheta*sphi*spsi   cpsi*sphi+ctheta*cphi*spsi  spsi*stheta;
            -spsi*cphi-ctheta*sphi*cpsi  -spsi*sphi+ctheta*cphi*cpsi cpsi*stheta;
            stheta*sphi                  -stheta*cphi                ctheta];        
   % Scaling matrix
   S = [ asq   0     0 ;
          0   bsq    0;
          0    0    csq;
       ];
 
   u= S^(-1)*A*(coord-repmat([x0 y0 z0].',1,size(coord,2)));
   tmp=rho*exp(-(u(1,:).^2+u(2,:).^2+u(3,:).^2)); 
   
  
   
   
   %tmp=rho*exp(-(u(1,:).^2+u(2,:).^2+u(3,:).^2)/2)/(sqrt(2*pi))^3 / det(S);   
   %tmp=rho*exp(-(u(1,:).^2+u(2,:).^2+u(3,:).^2)/2)/(sqrt(det(2*pi*S^2)));
   % To verify the normaliztions in the above formula make sure that 
   % sum(vol(:))*(2/(n-1))^3=1.

vol=vol+reshape(tmp,size(vol));
end

vol = reshape(vol,[n n n]);

return;


function out = flatten(in)

out = reshape(in,[1 numel(in)]);

return;
