function projections=cryo_project(volume,rot,n,precision)
%
% Project the given volume in a direction given by the rotations rot.
%
% Input parameters:
%   volume      3D array of size nxnxn of thevolume to project.
%   rot         Array of size 3-by-3-by-K containing the projection directions
%               of each projected image.
%   precision   Accuracy of the projection. 'single' or 'double'. Default
%               is single (faster).
%
% Output parameters:
%   projections 3D stack of size of projections. Each slice
%               projections(:,:,k) is a projection image of size nxn which
%               corresponds to the rotation rot(:,:,k). The number of
%               projections in the stack is equal to the number of
%               quaternions.
%
% The function supports creating projections whose size is different from
% the size of the proejcted volume. However, this code is still
% experimental and was not tested to verify that the resulting projections
% match the analytical ones to double precision.
%
% Note:
% To make the output of this function compatible with cryo_gen_projections
% call 
%   projections=permute(projections,[2 1 3]); 
% The FIRM reconstruction functions rely upon this permuted order. See
% 'examples' directory for detailed examples.
%
% Example:
%     voldef='C1_params';
%     rot = rand_rots(1);
%     n=129;
%     rmax=1;
%     vol=cryo_gaussian_phantom_3d(n,rmax,voldef);
%     p=cryo_project(vol,rot);
%     imagesc(p);
%
% Yoel Shkolnisky, August 2013.

if nargin<4
    precision='single';
end

if nargin<3
    n = size(volume,1);
end

if mod(n,2)==1
    range=-(n-1)/2:(n-1)/2;    
else
    range=-n/2+1/2:n/2-1/2;    
end

[I,J]=ndgrid(range,range);
I=I(:);
J=J(:);

N=numel(range); % If N is even, projection_fourier is not conjugate symmetric.
nv=size(volume,1);
% If the volume is too small for the given size, upsample it as needed.
if N>nv+1 % For compaibility with gen_projections, allow one pixel aliasing.
          % More precisely, it should be N>nv, however, by using nv+1 the
          % results match those of gen_projections.
    if mod(N-nv,2)==1
        error('Upsampling from odd to even sizes or vice versa is currently not supported');
    end
    dN=floor((N-nv)/2);
    fv=cfftn(volume);
    padded_volume=zeros(N,N,N,precision);
    padded_volume(dN+1:dN+nv,dN+1:dN+nv,dN+1:dN+nv)=fv;
    volume=icfftn(padded_volume);    
    assert(norm(imag(volume(:)))/norm(volume(:))<1.0e-5);
    nv=N; % The new volume size
end

K=size(rot, 3);
projections=zeros(N,N,K);

% Verify that we have only small imaginary components in the
% projcetions. How 'small' dependes on the accuracy. If the projections
% are 'single' allow for imaginary components that are about 10^-7. If
% 'double' allow for 10^-13.
if isa(volume,'single')
    imagtol=5.0e-7;
elseif isa(volume,'double')
    imagtol=5.0e-13;
else
    error('volume is not single nor double!?');
end


parfor k=1:K

%   if((mod(k,1000))==0)
% 		sprintf('%d projections done',k)
%   end
    R = rot(:,:,k);
    Rt=R.';
    
    % n_x, n_y, n_z are the image of the unit vectors in the x, y, z
    % directions under the inverse rotation
    
    n_x= Rt(:,1); %
    n_y= Rt(:,2);
    %n_z= Rt(:,3);  % not used - just for completeness
    
    
    P = I * n_x' + J * n_y';
    P= -2*pi*P/nv;
   
    projection_fourier = nufft3(volume, -P');
    
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

    if norm(imag(projection(:)))/norm(projection(:))>imagtol
        error('GCAR:imaginaryComponents','projection has imaginary components');
    end
    projection = real(projection);
    projections(:,:,k)=projection;
end

