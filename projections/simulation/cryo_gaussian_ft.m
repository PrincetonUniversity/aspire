function y=cryo_gaussian_ft(def,omega)
%
% Samples the analytic Fourier transform of the Gaussian phantom given by
% params at the points omega
%
% Yoel Shkolnisky, June 2007
%
% Y.S May 2013      Let the function take as a parameter the definition of
%                   the phantom.

if (~ismatrix(omega)) || (size(omega,2)~=3)
    error('omega must be an m by 3 array of frequnecies');
end

y=zeros(size(omega,1),1);
gparams=feval(def);
omega=omega.';

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
   
   u= S*A*omega;
   phase=exp(-1i*dot(omega,repmat([x0 y0 z0].',1,size(omega,2)),1));
   tmp=rho*phase.*exp(-(u(1,:).^2+u(2,:).^2+u(3,:).^2)./4).*(pi).^(3/2).*det(S);

   y=y+tmp.';
end
