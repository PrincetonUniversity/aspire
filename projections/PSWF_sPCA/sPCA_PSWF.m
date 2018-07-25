function [sPCA_data,denoised_images] = sPCA_PSWF_aspire(projections,nv,nfft_path,denoiseFlag,beta,T,useReflections)
% Compute steerable PCA on a dataset of images and perform denoising by
% truncation and shrinkage. The steps are as follows:
% 1) Expand images in PSWF basis
% 2) Compute SVD of each block of coefficients (for every angular index N)
% 3) The left singular vectors provide the steerable principal components
% 4) Detect componenets beyond the noise (by the Baik-Ben Arous-Peche phase transition)
% 4) (Optional) Shrink the expansion coefficients (in the basis of the steerable PC's) for improved de-noising (according to Gavish-Donoho (2017) -- note that singular-value shrinkage is equivalent to shrinkage of PCA expansion coefficients) 

% Input:    
%           projections: 3D array of images (third dimension enumerates over different images)
%           nfft_path: Path to NFFT by Potts
%           denoiseFlag: flag for applying shrinkage (improved denoising). Default = true.
%           beta: Oversampling factor (1 == no oversampleing, 0.5 == oversampling of x2, etc.). Default = 1. 
%           T: Truncation parameter (between 10^-6 and 10^6). With strong noise, use T ~ 10-10^3. Default = 10.
%           useReflections: Use reflections of images when computing SVD (compute singular values and vectors from real part of covariance only). Default = true.
%
% Output:   
%           sPCA_data: steerable PCA data structure

%% Default parameters
if ~exist('denoiseFlag','var')
    denoiseFlag=true;
end
if ~exist('beta','var')
    beta=1;
end
if ~exist('T','var')
    T=10;
end
if ~exist('useReflections','var')
    useReflections = true;
end
%% Add path
%addpath(genpath('./FPswfCoeffEval'));
%addpath(genpath(nfft_path));

%% Set params
remove_mean = 1;

%% Generate grids and sizes
sizeVec = size(projections); 
sizeX = sizeVec(1); sizeY = sizeVec(2); sizeZ = sizeVec(3);

nImages = sizeZ;

if mod(sizeX,2)==0
    evenOddFlag = 1;
else    % Make images always even-sized (consistentcy with previous SPCA version)  
%     evenOddFlag = 0;
    projections = projections(1:end-1,1:end-1,:);
    sizeVec = size(projections); 
    sizeX = sizeVec(1); sizeY = sizeVec(2); sizeZ = sizeVec(3);
    evenOddFlag = 1;
end

L=floor(sizeX/2);
if evenOddFlag==1
    [x, y]=meshgrid(-L:L-1, -L:L-1);    % Even-sized grid
    gridSize = 2*L;
else    
    [x, y]=meshgrid(-L:L,-L:L);     % Odd-sized grid
    gridSize = 2*L+1;
end

r=sqrt(x.^2+y.^2);
r_max = L;

%% Zero-pad if images are odd-sized
if evenOddFlag==0
    projections = [projections, zeros(gridSize,1,nImages)];
    projections = [projections; zeros(1,gridSize+1,nImages)];
    L = L+1;
    [x, y]=meshgrid(-L:L-1, -L:L-1);    % Even-sized grid
    r=sqrt(x.^2+y.^2);
    r_max = L;
end

%% Map images to PSWF expansion coefficients
[ PSWF_coeff, Psi, alpha, ang_freqs, ~, ~] = fastPswfCoeffEval( projections, x(r<=r_max)/r_max, y(r<=r_max)/r_max, r_max, beta, T, 64, 1 );

%% Perform Preliminary de-noising of PSWF coefficients by singular value shrinkage (for every angular index)
PSWF_coeff_denSVS = zeros(size(PSWF_coeff,1),size(PSWF_coeff,2));
rank = zeros(1,max(ang_freqs)+1);
wTot = [];
Psi_spca = [];
spca_coeff = [];
ang_freqs_spca = [];
rad_freqs_spca = [];
lambdaTot = [];
spca_coeff_denSVS = [];

if remove_mean
    mu = mean(PSWF_coeff(ang_freqs==0,:),2);
    PSWF_coeff(ang_freqs==0,:) = bsxfun(@minus,PSWF_coeff(ang_freqs==0,:),mu);
else
    mu = zeros(size(PSWF_coeff(ang_freqs==0,:),1),1);
end

for m=0:max(ang_freqs)
    % clc; 
    log_message(['Performing preliminairy de-noising of angular index: ',num2str(m),', out of ',num2str(max(ang_freqs))]);
    if m==0
        [PSWF_coeff_denSVS(ang_freqs==m,:), rank(m+1),w,pc,coeff,coeff_den,lambda] = matrixDenoise_sPCA(PSWF_coeff(ang_freqs==m,:),nv,false);
    else
        [PSWF_coeff_denSVS(ang_freqs==m,:), rank(m+1),w,pc,coeff,coeff_den,lambda] = matrixDenoise_sPCA(PSWF_coeff(ang_freqs==m,:),nv,useReflections);
    end
    wTot = [wTot w.'];
    Psi_spca = [Psi_spca Psi(:,ang_freqs==m)*pc];
    spca_coeff = [spca_coeff; coeff];
    ang_freqs_spca = [ang_freqs_spca; m*ones(rank(m+1),1)];
    rad_freqs_spca = [rad_freqs_spca; (1:rank(m+1)).'];
    lambdaTot = [lambdaTot lambda.'];
    spca_coeff_denSVS = [spca_coeff_denSVS; coeff_den];
end

%% Replace Prolates with steerable Principal Componenets
PSWF_N0 = Psi(:,ang_freqs==0);
Psi = Psi_spca;
ang_freqs = ang_freqs_spca;
if ~denoiseFlag    
    spca_coeff_denSVS = spca_coeff;
end

%% Form denoised images
if nargout>0
    I_denoise_svs = (PSWF_N0*mu)*ones(1,nImages) + Psi(:,ang_freqs==0)*spca_coeff_denSVS(ang_freqs==0,:) + 2*real(Psi(:,ang_freqs~=0)*spca_coeff_denSVS(ang_freqs~=0,:));
    denoised_images = zeros(2*L,2*L,nImages);
    for i=1:nImages
        tmpIm = zeros(2*L,2*L);
        tmpIm(r<=r_max) = real(I_denoise_svs(:,i));
        denoised_images(:,:,i) = tmpIm;
    end

%% Remove zero padding
    if evenOddFlag==0
        denoised_images = denoised_images(1:end-1,1:end-1,:);    
    end
end

%% Fill output structure
[sPCA_data.eigval,idx] = sort(lambdaTot,'descend');
sPCA_data.eigval = sPCA_data.eigval.';
sPCA_data.Freqs = ang_freqs(idx);
sPCA_data.RadFreqs = rad_freqs_spca(idx);
sPCA_data.Coeff = spca_coeff_denSVS(idx,:);
sPCA_data.Mean = mu;
sPCA_data.c = 0.5*beta;
sPCA_data.R = L;
sPCA_data.eig_im = zeros(4*L^2,numel(idx));
for i=1:numel(idx)
    tmpIm = zeros(2*L,2*L);
    tmpIm(r<=r_max) = Psi(:,idx(i));
    sPCA_data.eig_im(:,i) = tmpIm(:);
end
for i=1:size(PSWF_N0,2)
    tmpIm = zeros(2*L,2*L);
    tmpIm(r<=r_max) = PSWF_N0(:,i);
    sPCA_data.fn0(:,i) = tmpIm(:);
end
