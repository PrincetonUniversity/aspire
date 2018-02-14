function projections=cryo_project(volume,rot,n,precision,batch_size)
%
% Project the given volume in a direction given by the rotations rot.
%
% Input parameters:
%   volume      3D array of size nxnxn of thevolume to project.
%   rot         Array of size 3-by-3-by-K containing the projection directions
%               of each projected image.
%   precision   Accuracy of the projection. 'single' or 'double', or a
%               float respresenting the precision. Default is single
%               (faster).
%   batch_size  Process the images in batches of batch_size. Large batch
%               size is faster but requires more memory. Default is 100.
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
% Yoel Shkolnisky, February 2018.

if nargin<5
    batch_size=100;
end

if nargin<4
    precision=eps('single');
end

% precision can be 'single', 'double', or a float. The first two are to
% backward compatibility with other funtcions.
if ischar(precision) % So if we get 'single' or 'double', convert to float.
    precision=eps(precision);
end

% The function zeros gets as an argument 'single' or 'double'. We want to
% use 'single' when possible to save space. presicion_str is the
% appropriate string to pass to the functions 'zeros' later on.
if precision<eps('single')
    precision_str='single';
else
    precision_str='double';
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
    padded_volume=zeros(N,N,N,precision_str);
    padded_volume(dN+1:dN+nv,dN+1:dN+nv,dN+1:dN+nv)=fv;
    volume=icfftn(padded_volume);
    assert(norm(imag(volume(:)))/norm(volume(:))<1.0e-5);
    nv=N; % The new volume size
end

K=size(rot, 3);
batch_size=min(batch_size,K);
projection_batches=zeros(N,N,batch_size,ceil(K/batch_size),precision_str);
% Each batch of images is processed and stored independently. All batches
% will be merged below into a single output array.

% Verify that we have only small imaginary components in the
% projcetions. How 'small' dependes on the accuracy.
imagtol=precision*5;

% Accuracy parameter for NUFFT below.
nufft_opt.epsilon=precision;

parfor batch=1:ceil(K/batch_size)
    % It may be that the number of remained images is less than batch_size.
    % So compute the actual_batch_size.
    actual_batch_size=min(batch_size,K-(batch-1)*batch_size);
    P=zeros(numel(I)*actual_batch_size,3); % Sampling points in Fourier
    % domain for all images of the current batch.
    startidx=(batch-1)*batch_size; % The linear index of the first image in
    % the current batch of images.
    
    for k=1:actual_batch_size
        %   if((mod(k,1000))==0)
        % 		sprintf('%d projections done',k)
        %   end
        R = rot(:,:,startidx+k);
        Rt=R.';
        
        % n_x, n_y, n_z are the image of the unit vectors in the x, y, z
        % directions under the inverse rotation
        
        n_x= Rt(:,1); %
        n_y= Rt(:,2);
        %n_z= Rt(:,3);  % not used - just for completeness
        
        
        P((k-1)*numel(I)+1:k*numel(I),:) = I * n_x' + J * n_y';
    end
    P= -2*pi*P/nv;
    
    % NUFFT all images in the current batch
    projection_fourier = nufft3(volume, -P',nufft_opt);
    projection_fourier=reshape(projection_fourier,numel(I),actual_batch_size);
    P = reshape(P, [numel(I), actual_batch_size, 3]);
    
    if mod(n,2)==0
        projection_fourier = projection_fourier.*exp(1i.*sum(P,3)./2);
        Irep=repmat(I,1,actual_batch_size);
        Jrep=repmat(J,1,actual_batch_size);
        projection_fourier = projection_fourier ...
            .* exp(2*pi*1i*(Irep+Jrep-1)/(2*n));
    end
    
    % IFFT the images from Fourier space to real space.
    projections_temp=zeros(N,N,batch_size,precision_str);
    for k=1:actual_batch_size
        temp = reshape(projection_fourier(:,k), N, N);
        temp = ifftshift(temp);
        projection = fftshift(ifft2(temp));
        
        if mod(n,2)==0
            projection = projection .* reshape(exp(2*pi*1i*(I+J)/(2*n)),n,n);
        end
        
        if norm(imag(projection(:)))/norm(projection(:))>imagtol
            error('GCAR:imaginaryComponents','projection has imaginary components');
        end
        projection = real(projection);
        projections_temp(:,:,k)=projection;
    end
    projection_batches(:,:,:,batch)=projections_temp;
end

% Merge projection_batches into a single output array
projections=zeros(N,N,K,precision_str);
for batch=1:ceil(K/batch_size)
    actual_batch_size=min(batch_size,K-(batch-1)*batch_size);
    startidx=(batch-1)*batch_size;
    projections(:,:,startidx+1:startidx+actual_batch_size)=...
        projection_batches(:,:,1:actual_batch_size,batch);
end
end
