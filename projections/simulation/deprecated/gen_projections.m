function [projections, noisy_projections, shifts, q] = ...
    gen_projections(K,SNR,max_shift,step_size,noise_type,silent,fprecomp)
%
% Generate K projections with the given SNR and radnom shifts.
%
% SNR is the signal-to-noise-ratio of each projection, defined as the
% variance of the clean projection divded by the variance of the noise.
% The default SNR is 1/60.
%
% max_shift is the maximal random shift (in pixels) introduced to each
% projection. max_shift must be an integer and the resulting random shifts
% are integers between -max_shift to +max_shift. The default is 0 (no
% shift). 
%
% step_size is the resolution used to generate shifts. step_size=1 allows
% for all integer shifts from -max_shift to max_shift (in both x and y 
% directions). step_size=2 generates shifts between -max_shift and
% max_shift with steps of 2 pixels, that is, shifts of 0 pixels, 2 pixels,
% 4 pixels, and so on. Default step_size is 3.
%
% q are the quaternions used to generate the projections.
%
% noise_type can be either 'gaussian' or 'colored'.
%
% Set silent to 1 to eliminate progess messages.
%
% fprecomp is a filename that contains precomputed data. Use this option
% whenever you need the same projections with different noise levels. 
% The clean shifted projections are loaded from the file and noise is
% added. If precomp is specified but does not exist, the generated
% projections are saved into the specified file.
%
% This function uses accurate interpolation based on non-equally spaced FFT
% to generete the projections. gen_projections_v3 a crude interpolation
% which resulted in larger discrepancies between common-lines.
%
% Shifts is a Kx2 array of the 2D shift of each projection. shifts(k,1)=1
% means that the k projection was shifted to the right by 1 pixel.
% shifts(k,2)=1 means that the projection was shifted down by 1 pixel.
%
% Revisions:
%   02/03/2009   Filename changed from gen_projections_v4.m to
%                gen_projections.m
%   10/19/2010   Use parfor for speedup.
%
% Yoel Shkolnisky and Amit Singer, January 2009.

initstate;
%rng('default')
%rng shuffle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate K projections of a given volume
% L is the number of radial lines in each projection
% The projections are random
% In order to quickly generate the projections
% we use the Fourier-Projection Slice Theorem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    SNR=1/60;
    fprintf('Using default SNR=%4.3f\n',SNR);
end

if nargin<3
    max_shift=0; % do not shift
    fprintf('Using default max_shift=%d\n',max_shift);
end

if nargin<4
    step_size=3;
    fprintf('Using default step_size=%d\n',step_size);
end

if nargin<5
    noise_type='gaussian';
    fprintf('Using Gaussian noise\n');
end

if nargin<6
    silent=0;
end

if step_size<0
    warning('GCAR:invalidArgument','step_size cannot be zero or negative. Setting step_size to 0...');
    step_size=0;
end

%% Load precomp data if needed and verify it against the given parameters

precomp=0;
if exist('fprecomp','var')       
    [pathstr, name, ext] = fileparts(fprecomp);
    if isempty(ext)
        fprecomp=[fprecomp,'.mat'];
    end
    if exist(fprecomp,'file')
        if ~silent
            fprintf('Loading precomputed data...');
        end
        projections=0; % To eliminate warning when referring to this variable.
        load(fprecomp); % loads projections, noisy_projections, shifts, q,and N.
        precomp=1;
        if ~silent
            fprintf('Finshed!\n');
        end

        if ~exist('projections','var') || ~exist('shifts','var') ||...
                ~exist('q','var') || ~exist('N','var')
            warning('GCAR:incompatibleArgument',...
                'Saved data is corrupted. Computing projections...');
            precomp=0;
        elseif size(projections,3)~=K
            warning('GCAR:incompatibleArgument',...
                'Saved data does not much given parameters. Computing projections...');
            precomp=0;
        end
    end
end
if (precomp==0)
    
    if ~silent
        fprintf('Loading volume...');
    end
    load cleanrib.mat; % This is the clean E. Coli Ribosome

% Another volume
%    load vol129

% C2 volume
      %load ~/data/work/molec/datasets/AMPA/C2CenteredVol.mat
      %load ~/data/work/molec/datasets/AMPA/Map3A
      %volref=Map3A;
      %volref=volref(1:63,1:63,1:63);
      
%      load reconstructed_vol_clean_K200C2rib.mat
%      volref=volr;
%      volref=GaussFilt3(volref,0.1);
     
    if ~silent
        fprintf('Finshed!\n');
        fprintf('Generating Projections...');
    end
    
    %volref=permute(volref,[2 1 3]);
    volume = volref;

    n = size(volume,1);

    %%%%%% Generate K random quaternions (a,b,c,d) %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                                                                     %%
    % The 3-sphere S^3 in R^4 is a double cover of the rotation group SO(3)
    % SO(3) = RP^3
    % We identify unit norm quaternions a^2+b^2+c^2+d^2=1 with group elements
    % The antipodal points (-a,-b,-c,-d) and (a,b,c,d) are identified as the
    % same group elements, so we take a>=0.

    q = randn(4,K);

    % q(:,1) = [0;1;0;1];

    l2_norm = sqrt(q(1,:).^2 + q(2,:).^2 + q(3,:).^2 + q(4,:).^2);

    for i=1:4;
        q(i,:) = q(i,:) ./ l2_norm;
    end;

    for k=1:K;
        if (q(1,k) < 0)
            q(:,k) = -q(:,k);
        end;
    end;

    %%%%%% Convert quaternions into 3 x 3 rotation matrices %%%%%%%%%%%%%%%%%
    %                                                                     %%
    % Using the Euler formula
    rot_matrices = zeros(3,3,K);

    for k=1:K;
        rot_matrix = q_to_rot(q(:,k));
        rot_matrices(:,:,k) = rot_matrix;
    end;

    % calculate inverse rotation matrices (just transpose)
    inv_rot_matrices = zeros(3,3,K);

    %inv_rot_matrix = zeros(3);
    for k=1:K;
        rot_matrix = rot_matrices(:,:,k);
        inv_rot_matrix = rot_matrix'; % inv(R)=R^T
        inv_rot_matrices(:,:,k) = inv_rot_matrix;
    end;

    %%%%%% Calculate great circles over S^2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % The circle is uniquely described by the orthogonal system (n_x,n_y,n_z)
    % n_z is perpendicular to the plane of the circle,
    % which is given by the two vectors n_x and n_y
    %
    % n_x, n_y, n_z are the image of the unit vectors in the x, y, z
    % directions under the inverse rotation

    n_x(:,:) = inv_rot_matrices(:,1,:); % n_x = zeros(3,K);
    n_y(:,:) = inv_rot_matrices(:,2,:);
    %n_z(:,:) = inv_rot_matrices(:,3,:);  % not used - just for completeness

    len=-fix(n/2):fix(n/2);
    N=numel(len); % If N is even, projection_fourier is not conjugate symmetric.
    [I,J]=meshgrid(len,len);
    I=I(:);
    J=J(:);
%    projections = zeros(2*N-1,2*N-1,K);
    projections = zeros(N,N,K);

    % Go concurrent
    poolreopen; 

    parfor k=1:K;
        P = I * n_x(:,k)' + J * n_y(:,k)';
        P= -2*pi*P/(n+1);
        
        projection_fourier = nufft3(volume, -P);
        projection_fourier = reshape(projection_fourier, N, N);

        % zero pad for upsampling
%        padded_projection=zeros(2*N-1);
%        padded_projection((N-1)/2+1:3*(N-1)/2+1,(N-1)/2+1:3*(N-1)/2+1)=projection_fourier;
%        projection_fourier=padded_projection;

        projection_fourier = ifftshift(projection_fourier);
        projection = fftshift(ifft2(projection_fourier));

        if norm(imag(projection(:)))/norm(projection(:)) >1.0e-8
            warning('GCAR:imaginaryComponents','projection has imaginary components');
        end

        projections(:,:,k) = projection;
    end;

    projections = real(projections);


    % shift projections
    shifts=zeros(K,2);
    if step_size>0
        shifts=round((rand(K,2)-1/2)*2*max_shift/step_size)*step_size;
    end

    if exist('fprecomp','var')
        save(fprecomp,'projections','shifts','q','N','len','n');
    end
    
    if ~silent
        fprintf('Finished!\n');
    end
    
end

% Shift projections

if ~silent
    fprintf('Adding shifts...');
end

%[omega_x,omega_y]=ndgrid(-(N-1):N-1,-(N-1):N-1);
%omega_x=-2*pi.*omega_x/(2*N-1); omega_y=-2*pi.*omega_y/(2*N-1);
[omega_x,omega_y]=ndgrid(len,len);
omega_x=-2*pi.*omega_x/(n+1); omega_y=-2*pi.*omega_y/(n+1);


parfor k=1:K
    p=projections(:,:,k);
    pf=fftshift(fft2(ifftshift(p)));
    phase_x=omega_x.*shifts(k,1);
    phase_y=omega_y.*shifts(k,2);
    pf=pf.*exp(sqrt(-1)*(phase_x+phase_y));
    p2=fftshift(ifft2(ifftshift(pf)));
    projections(:,:,k)=real(p2);
end

if ~silent
     fprintf('Finished!\n');
end


if ~silent
    fprintf('Adding noise...');
end

% Add noise to projections
p = size(projections, 1);
noisy_projections=zeros(size(projections));

rand('state',1234);
randn('state',1234);

if ~strcmpi(noise_type,'gaussian')
    % If we need colored noise, the load the appropriate filter.
    load color_filter.mat;
    noise_response=ifftshift(S);  %for optimization, so we can use fft2 and
%ifft2 below instead of cfft2 and icfft2
end

lowidx=-(p-1)/2+p+1;
highidx=(p-1)/2+p+1;

for k=1:K
    proj=projections(:,:,k);
    sigma=sqrt(var(proj(:))/SNR);
    gn=randn(2*p+1);

    if strcmpi(noise_type,'gaussian')
        cn=gn;
    else
        cn=ifft2(fft2(gn).*noise_response);
    end
    cn=cn(lowidx:highidx,lowidx:highidx);
    cn=cn/std(cn(:));
    cn=cn.*sigma;
    noisy_projections(:,:,k) = proj + cn;
end

% for k=1:K
%     proj=projections(:,:,k);
%     sigma=sqrt(var(proj(:))/SNR);
%     noise=randn(p,p) * sigma;
%     noisy_projections(:,:,k) = proj + noise;
% end

if ~silent
    fprintf('Finished!\n');
end
