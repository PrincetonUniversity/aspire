function p=cryo_project_gaussian(def,n,rmax,rots)
%
% Compute 2D analytic projections of the 3D Guassian phantom.
%
% rots is an array of rotation matrices, where each matrix corresponds to a
% projection direction. Each projection is of size nxn, sampled on a
% Cartesian grid between -rmax and rmax (default [-1,1]). 
%
% Input parameters:
%     def   Filename of parameters file of the phantom.
%     n     Dimension of each projection.
%     rmax  Physical dimensions of each projection. Each projection has n
%           samples in each dimension between -rmax and rmax.
%     rots  Array of size 3-by-3-by-K containing rotation matrices along which
%           to project.
%
%
% Output parameters;
%     p         3D array of size n x n x size(rots, 3).
%               p(:,:,k) is the projection in the direction determined by
%               the rotation matrix rot(:,:,k).
%
% Example:
%
%     rots = rand_rots(100);
%     p=cryo_project_gaussian('C1_params',65,1,rots);
%
% Another example: (code sanity check)
%     rot = eye(3);
%     n=129; 
%     rmax=1;
%     vol=cryo_gaussian_phantom_3d('C4_params',n,rmax); 
%     p=cryo_project_gaussian('C4_params',n,rmax,rot);
%     p2=sum(vol,3)*2/(n-1);
%     max(abs(p2(:)-p(:)))  % should be tiny
%     norm(p2(:)-p(:))/norm(p(:))  % should be tiny
%
% Yoe Shkolnisky, August 2013
% Based on gaussian_projections3 from June 2007
%

if nargin<4
    rmax=1;
end

%t=linspace(-abs(rmax),abs(rmax),n);
if mod(n,2)==1
    t=(-(n-1)/2:(n-1)/2) / ((n-1)/2)*rmax; 
else
    t=(-n/2+1/2:n/2-1/2)/(n/2)*rmax;
end

[u,v]=ndgrid(t,t); % Do not use ndgrid to be compatible with gen_projections.

s_u=size(u);

u=u(:);
v=v(:);
l_u=size(u,1);

nproj=size(rots,3);
p=zeros(s_u(1),s_u(2),nproj);
gparams = feval(def);

for j = 1:size(gparams,1)
    rho = gparams(j,1);          % Amplitude change for this Gaussian
    asq = gparams(j,2);          % std x
    bsq = gparams(j,3);          % std y
    csq = gparams(j,4);          % std z
    x0 = gparams(j,5);           % x offset
    y0 = gparams(j,6);           % y offset
    z0 = gparams(j,7);           % z offset
    r_phi = gparams(j,8)*pi/180;   % first Euler angle in radians
    r_theta = gparams(j,9)*pi/180; % second Euler angle in radians
    r_psi = gparams(j,10)*pi/180;  % third Euler angle in radians

    cphi = cos(r_phi);
    sphi = sin(r_phi);
    ctheta = cos(r_theta);
    stheta = sin(r_theta);
    cpsi = cos(r_psi);
    spsi = sin(r_psi);

    % Euler rotation matrix
    A = [cpsi*cphi-ctheta*sphi*spsi   cpsi*sphi+ctheta*cphi*spsi  spsi*stheta;
        -spsi*cphi-ctheta*sphi*cpsi  -spsi*sphi+ctheta*cphi*cpsi cpsi*stheta;
        stheta*sphi                  -stheta*cphi                ctheta];

    % Scaling matrix
    S = [asq   0     0;
          0   bsq    0;
          0    0    csq;
        ];

    r0=[x0;y0;z0];
    
    inv_rot_matrices = permute(rots, [2 1 3]);
    
    for k=1:nproj
        %
        % orthogonal basis for the projection direction.
        % tau is the projection direction. alpha and beta are the x-axis
        % and y-axis in the plane perpendicular to the projection
        % direction. The triplet [alpha,beta,tau] is a right-handed system.
        %
%         tau=[cos(theta(k))*cos(phi(k)); sin(theta(k))*cos(phi(k)); sin(phi(k))];
%         alpha=[sin(theta(k)); -cos(theta(k)); 0];
%         beta=[cos(theta(k))*sin(phi(k)); sin(theta(k))*sin(phi(k)); -cos(phi(k))];
        
        alpha = inv_rot_matrices(:,1,k); 
        beta  = inv_rot_matrices(:,2,k);
        tau   = inv_rot_matrices(:,3,k);

        % The system should be right-handed
        if abs(det([alpha beta tau])-1)>1.0e-6
            warning('System is not right-handed');
        end

        % Variables that end with _t refer to a coordinate systems that
        % has underwent the transformation S^(-1)*A

        tau_t=S^(-1)*A*tau;
        gamma=norm(tau_t);
        tau_t=tau_t/gamma;

        phi_t=asin(tau_t(3));

        if abs(abs(tau_t(3))-1)<1.0e-11
            theta_t=pi/2;
        else
            theta_t=atan2(tau_t(2),tau_t(1));
        end

        alpha_t=[sin(theta_t); -cos(theta_t); 0];
        beta_t=[cos(theta_t)*sin(phi_t); sin(theta_t)*sin(phi_t); -cos(phi_t)];

        w=(u*(alpha.')+v*(beta.')).'-repmat(r0,1,l_u);
        w=S^(-1)*A*w;        
        u0=alpha_t.'*w;
        v0=beta_t.'*w;

        %tmp=(rho/gamma)/(2*pi)/det(S)*exp(-(u0.^2+v0.^2)/2);
        tmp=(rho/gamma)*sqrt(pi)*exp(-(u0.^2+v0.^2));
        p(:,:,k)=p(:,:,k)+reshape(tmp,s_u);

    end
end


