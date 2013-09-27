function projections=cryo_project(volume,q,n,precision)
%
% Project the given volume in a direction given by the quaternions q.
%
% Input parameters:
%   volume      3D array of size nxnxn of thevolume to project.
%   q           List of quaternions determining the projection direction of
%               each projected image.
%   precision   Accuracy of the projection. 'single' or 'double'. Default
%               is double (slower).
%
% Output parameters:
%   projections 3D stack of size of projections. Each slice
%               projections(:,:,k) is a projection image of size nxn which
%               corresponds to the quaternion q(:,k). The number of
%               projections in the stack is equal to the number of
%               quaternions.
%
% Example:
%     voldef='C1_params';
%     q = qrand(1);
%     n=129;
%     rmax=1;
%     vol=cryo_gaussian_phantom_3d(n,rmax,voldef);
%     p=cryo_project(vol,q);
%     imagesc(p);
%
% Yoel Shkolnisky, August 2013.

if nargin<4
    precision='double';
end

if nargin<3
    n = size(volume,1);
end

if mod(n,2)==1
    range=-(n-1)/2:(n-1)/2;    
else
    range=-n/2+1/2:n/2-1/2;    
end

N=numel(range); % If N is even, projection_fourier is not conjugate symmetric.
nv=size(volume,1);
[I,J]=ndgrid(range,range);
I=I(:);
J=J(:);
prepdata=nufft_t_3d_prepare_2(volume,precision);

K=size(q,2);
projections=zeros(N,N,K);

poolreopen;
parfor k=1:K
    
    R=q_to_rot(q(:,k));
    Rt=R.';
    
    % n_x, n_y, n_z are the image of the unit vectors in the x, y, z
    % directions under the inverse rotation
    
    n_x= Rt(:,1); %
    n_y= Rt(:,2);
    %n_z= Rt(:,3);  % not used - just for completeness
    
    
    P = I * n_x' + J * n_y';
    P= -2*pi*P/nv;
   
    
    %    projection_fourier = nufft_t_3d(volume,P,'single');
    projection_fourier = nufft_t_3d_execute_2(P,prepdata);
    
    if mod(n,2)==0
        projection_fourier = projection_fourier.*exp(1i.*sum(P,2)./2);
        projection_fourier = projection_fourier .* exp(2*pi*1i*(I+J-1)/(2*n));
    end
    projection_fourier = reshape(projection_fourier, N, N);       
    projection_fourier = ifftshift(projection_fourier);
    projection = fftshift(ifft2(projection_fourier));

    if mod(n,2)==0
        projection = projection .* reshape(exp(2*pi*1i*(I+J)/(2*n)),n,n);
    end

    
    if norm(imag(projection(:)))/norm(projection(:)) >1.0e-8
        error('GCAR:imaginaryComponents','projection has imaginary components');
    end
    projection = real(projection);
    projections(:,:,k)=projection;
end