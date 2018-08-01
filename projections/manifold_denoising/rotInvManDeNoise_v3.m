function [ projections_man_den,projections_svs_den,projections_sPca,projections_man_EM_den,rndPermIdx,lScore] = rotInvManDeNoise_v3( filePath,beta,T,N,nImSMH,maxEigIdx,useSPCA,useWienerFilt,shuffle,nImCV,olN,nn_p)
%
% XXXFUNCTIONAMEXXX     Manifold denoising of cryo-EM images.
%
% Input arguments:
%   filePath     Filename in MRC format of noisy projections.
%   beta         PSWFs oversampling ratio factor. beta=1 mean the images
%         were sampled at Nyquist, beta=0.5 means they were sampled at 1/2
%         Nyquist and so on. Default is beta=1. 
%   T            PSWFs truncation parameter. For definition see Applied and
%         Computational Harmonic Analysis 43(3):381-403, November 2017.
%         Default is T=10. 
%   N            Number of images to use from the dataset. Default is all
%         images.
%   nImSMH       Number of images to be used in the steerable graph
%         Laplacian. Default nImSMH=5000. 
%   maxEigIdx    Maximal number of eigenfunctions to take from the
%         steerable graph Laplacian. Default maxEigIdx=2000. 
%   useSPCA      Replace PSWF basis and coefficients with sPCA basis and
%         coefficients. Default useSPCA=1.
%   useWienerFilt  Use Wiener filtering when computing pair-wise distances 
%         in the steerable graph Laplacian. Default useWienerFilt=1.
%   shuffle      Randomely shuffle the images. Default shuffle=1.
%   nImCV        Number of images to be used in cross validation to find
%                optimal parameters (half of nImCV is used for train i.e.
%                denoising, and the other half for likelihood test).
%                Default nImCV=10000.
%   olN          Number of images to classify as outliers and remove (based
%         on pixel variance extermum). Default is to remove 25% of the
%         input images.
%   nn_p        Percentage of nearest-neighbors to retain when constructing
%         the steerable graph Laplacian. Default nn_p=5%.
%
%   You can omit any optional input argument or pass an empty array []. In
%   such a case a default value will be used.
%
% Output arguments:
%   projections_man_den     Denoised images using manifold denosing.
%   projections_svs_den     Denoised images using steerable-PCA and
%         shrinkage.
%   projections_sPca        Denoised images using steerable PCA without
%         shrinkage.
%   projections_man_EM_den  Denoised images using manifold denoising
%         where XXX                
%   rndPermIdx              Permuation applied to the input images
%   lScore                  XXX
%
% For example see XXX.m
%
% Boris Lande, May 2018.

if ~exist('beta','var') || isempty(beta)
    beta=1;
end

if ~exist('T','var') || isempty(T)
    T=10;
end

if ~exist('N','var') || isempty('N')
    stack=imagestackReader(filePath); 
    N=stack.dim(3);
end

if ~exist('nImSMH','var') || isempty(nImSMH)
    nImSMH=5000;
end

if ~exist('maxEigIdx','var') || isempty(maxEigIdx)
    maxEigIdx=2000;
end

if ~exist('useSPCA','var') || isempty(useSPCA)
    useSPCA=1;
end

if ~exist('useWienerFilt','var') || isempty(useWienerFilt)
    useWienerFilt=1;
end

if ~exist('shuffle','var') || isempty(shuffle)
    shuffle=1;
end

if ~exist('nImCV','var') || isempty(nImCV)
    nImCV=10000;
end

if ~exist('olN','var') || isempty(olN)
    olN=round(N*0.1/4);
end

if ~exist('nn_p','var') || isempty(nn_p)
    nn_p=5;
end

if ~is_gpu
    error('CUDA enabled GPU is required');
end


%% Generate grids and sizes
stack=imagestackReader(filePath); 
sizeVec = stack.dim; 
sizeX = sizeVec(1); sizeY = sizeVec(2); sizeZ = sizeVec(3);

if ~isempty(N)
    nImages = min(N,sizeZ);
else
    nImages = 2 * sizeZ;
end
nImSMH = min(nImSMH,nImages);

L=floor(min(sizeX,sizeY)/2);
[x, y]=meshgrid(-L:L-1, -L:L-1);    % Even-numbered grid
% [x, y]=meshgrid(-L:L,-L:L); 

r=sqrt(x.^2+y.^2);
r_max = L;

%% Print parameters
log_message('Starting manifold denoising on %s',filePath);
log_message('Denoising parameters:');
log_message('\t beta=%4.2f',beta);
log_message('\t T=%d',T);
log_message('\t N=%d',nImages);
log_message('\t nImSMH=%d',nImSMH);
log_message('\t maxEigIdx=%d',maxEigIdx);
log_message('\t useSPCA=%d',useSPCA);
log_message('\t useWienerFilt=%d',useWienerFilt);
log_message('\t shuffle=%d',shuffle);
log_message('\t nImCV=%d',nImCV);
log_message('\t olN=%d',olN);
log_message('\t nn_p=%d',nn_p);
gpudev=gpuDevice;
log_message('\t gpuDeviceIdx=%d',gpudev.Index);

%% Load data
sVarVec = zeros(1,nImages);
nVarVec = zeros(1,nImages);
projections = zeros(2*L,2*L,nImages);
log_message('Loading data...');
for j = 1:nImages
        currImage=stack.getImage(j);
        currImage = currImage(1:2*L,1:2*L);
%         currImageT = currImage.';
%         currImage = currImage(1:2*L,1:2*L)/sqrt(var(currImage(r>r_max)));
        projections(:,:,j) = currImage;
%         projections(:,:,j + nImages / 2) = currImageT;
        sVarVec(j) = var(currImage(r<=r_max));
%         sVarVec(j + nImages / 2) = var(currImageT(r<=r_max));
        nVarVec(j) = var(currImage(r>r_max));
%         nVarVec(j + nImages / 2) = var(currImageT(r>r_max));
end

%% Remove outliers
[~,idx] = sort(nVarVec,'ascend');
idx_ol1 = [idx(1:olN),idx((end-olN+1):end)];
projections(:,:,idx_ol1) = [];
nVarVec(idx_ol1) = [];
sVarVec(idx_ol1) = [];

[~,idx] = sort(sVarVec,'ascend');
idx_ol2 = [idx(1:olN),idx((end-olN+1):end)];
projections(:,:,idx_ol2) = [];
nVarVec(idx_ol2) = [];
sVarVec(idx_ol2) = [];

nImages = size(projections,3);

%% Shuffle the images
if shuffle
    rndPermIdx = randperm(nImages);
    projections = projections(:,:,rndPermIdx);
    nVarVec = nVarVec(rndPermIdx);
    sVarVec = sVarVec(rndPermIdx);
else
    rndPermIdx = 1:nImages;
end

%% Estimate noise variance and signal power
nv = mean(nVarVec);
sp = mean(sVarVec);

snr = 10*log10(sp/nv-1);
log_message('Estimated data SNR %5.1fdB (%4.2e)',snr,10^(snr/10));

%% Map images to PSWF expansion coefficients
log_message('Computing PSWF expansion coefficients...');
[ PSWF_coeff, Psi, ~, ang_freqs, ~, ~] = fastPswfCoeffEval( projections, x(r<=r_max)/r_max, y(r<=r_max)/r_max, r_max, beta, T, 64, 1 );

%% Perform Preliminary de-noising of PSWF coefficients by singular value shrinkage (for every angular index)
log_message('Computing steerable PCA and performing preliminary denoising...')
PSWF_coeff_denSVS = zeros(size(PSWF_coeff,1),size(PSWF_coeff,2));
rank = zeros(1,max(ang_freqs)+1);
wTot = [];
Psi_spca = [];
spca_coeff = [];
ang_freqs_spca = [];
sdTot = [];

mu =mean(PSWF_coeff(ang_freqs==0,:),2);
PSWF_coeff(ang_freqs==0,:) = bsxfun(@minus,PSWF_coeff(ang_freqs==0,:),mu);
mu = Psi(:,ang_freqs==0)*mu;

for m=0:max(ang_freqs)
%     disp(['Performing preliminairy de-noising of angular index: ',num2str(m),', out of ',num2str(max(ang_freqs))]);
    [PSWF_coeff_denSVS(ang_freqs==m,:), rank(m+1),w,pc,coeff,sd_c] = matrixDenoise(PSWF_coeff(ang_freqs==m,:),nv);
    wTot = [wTot w.'];
    Psi_spca = [Psi_spca Psi(:,ang_freqs==m)*pc];
    spca_coeff = [spca_coeff; coeff];
    ang_freqs_spca = [ang_freqs_spca; m*ones(rank(m+1),1)];
    sdTot = [sdTot sd_c.'];
end

%% Replace Prolates with steerable Principal Componenets
if useSPCA
    PSWF_coeff = spca_coeff;    
    Psi = Psi_spca;
    ang_freqs = ang_freqs_spca;
    PSWF_coeff_denSVS = diag(wTot)*spca_coeff;
end

%% View preliminary de-noising results for debugging
% testIdx = 1;
% I_denoise_pre = mu + Psi(:,ang_freqs==0)*PSWF_coeff_denSVS(ang_freqs==0,:) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff_denSVS(ang_freqs~=0,:));
% image_d_disp = zeros(2*L,2*L);
% image_d_disp(r<=r_max) = real(I_denoise_pre);
% figure; imagesc(image_d_disp); colormap(gray)

% I_denoise_svs = mu*ones(1,nImages) + Psi(:,ang_freqs==0)*PSWF_coeff_denSVS(ang_freqs==0,:) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff_denSVS(ang_freqs~=0,:));
% projections_svs_den = zeros(2*L,2*L,nImages);
% for i=1:nImages
%     tmpIm = zeros(2*L,2*L);
%     tmpIm(r<=r_max) = real(I_denoise_svs(:,i));
%     projections_svs_den(:,:,i) = tmpIm;
% end
% figure; viewstack(projections_svs_den(:,:,1:1000),10,10,0);

%% Do mini-cross-validation on a subset of the data (over parameters lambda_c and epsilon)
log_message('Starting cross-validation procedure for optimal parameters...')
% nTheta = 256;
nTheta = 512;
nImCV = min(nImCV,nImages);
trainSize = round(nImCV/2);
nn_cv=ceil(nn_p*trainSize/100);

if useWienerFilt
    x_sMH = PSWF_coeff_denSVS(:,1:trainSize);
    x = PSWF_coeff(:,1:nImCV);
    W_eps = sqrt(sum(wTot.^4))*nv;
else
    x_sMH = PSWF_coeff(:,1:trainSize);
    x = PSWF_coeff(:,1:nImCV);
    W_eps = 1*sqrt(sum(rank))*nv;
end

lambda_c = 0.5;

empLkh = [];
cnt = 1;
% [vCell,dCell,~] = evalSMH(x_sMH,ang_freqs,W_eps,nTheta,[],t_Thr,0);
[vCell,dCell] = evalSMH_sparse(x_sMH,ang_freqs,W_eps,nTheta,maxEigIdx,nn_cv,0);
empLkh(cnt) = getCVlogL_gpu(x,ang_freqs,nTheta,trainSize,nv,vCell,dCell,lambda_c);

%disp(['Iteration ',num2str(cnt),': lambda_c = ',num2str(lambda_c),', epsilon = ',num2str(W_eps), ', Log-likelihood = ',num2str(empLkh(cnt))]);
log_message('Iteration %d: lambda_c = %4.2f, epsilon = %4.2f, Log-likelihood = %4.2f',cnt,lambda_c,W_eps,empLkh(cnt));
dlMCV = 0.05;
dEps = sqrt(2);
endFlag = 0;
while(1)        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Do lambda_c coordinate iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % - Determine direcion 
    empLkh_inc = getCVlogL_gpu(x,ang_freqs,nTheta,trainSize,nv,vCell,dCell,lambda_c+dlMCV); % Determine direcion
    if empLkh_inc>=empLkh(cnt)
        cnt = cnt + 1;
        empLkh(cnt) = empLkh_inc;
        lambda_c = lambda_c + dlMCV;
        %disp(['Iteration ',num2str(cnt),': lambda_c = ',num2str(lambda_c),', epsilon = ',num2str(W_eps), ', Log-likelihood = ',num2str(empLkh(cnt))]);
        log_message('Iteration %d: lambda_c = %4.2f, epsilon = %4.2f, Log-likelihood = %4.2f',cnt,lambda_c,W_eps,empLkh(cnt));        
    else
        dlMCV = -dlMCV; % Change direction
    end    
    % - Run in that direcion 
    while(1)
        empLkh_inc = getCVlogL_gpu(x,ang_freqs,nTheta,trainSize,nv,vCell,dCell,lambda_c+dlMCV);
        if empLkh_inc>=empLkh(cnt)
            cnt = cnt + 1;
            empLkh(cnt) = empLkh_inc;
            lambda_c = lambda_c + dlMCV;
            %disp(['Iteration ',num2str(cnt),': lambda_c = ',num2str(lambda_c),', epsilon = ',num2str(W_eps), ', Log-likelihood = ',num2str(empLkh(cnt))]);
            log_message('Iteration %d: lambda_c = %4.2f, epsilon = %4.2f, Log-likelihood = %4.2f',cnt,lambda_c,W_eps,empLkh(cnt));
        else
            break;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Do W_eps coordinate iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    first_W_eps_update = 1;
    % - Determine direcion 
%     [vCell_inc,dCell_inc,~] = evalSMH(x_sMH,ang_freqs,dEps*W_eps,nTheta,[],t_Thr,0);  
    [vCell_inc,dCell_inc] = evalSMH_sparse(x_sMH,ang_freqs,dEps*W_eps,nTheta,maxEigIdx,nn_cv,0);
    empLkh_inc = getCVlogL_gpu(x,ang_freqs,nTheta,trainSize,nv,vCell_inc,dCell_inc,lambda_c);     
    if empLkh_inc>=empLkh(cnt)
        first_W_eps_update = 0;
        cnt = cnt + 1;
        empLkh(cnt) = empLkh_inc;
        W_eps = dEps*W_eps;
        %disp(['Iteration ',num2str(cnt),': lambda_c = ',num2str(lambda_c),', epsilon = ',num2str(W_eps), ', Log-likelihood = ',num2str(empLkh(cnt))]);
        log_message('Iteration %d: lambda_c = %4.2f, epsilon = %4.2f, Log-likelihood = %4.2f',cnt,lambda_c,W_eps,empLkh(cnt));
        vCell = vCell_inc;
        dCell = dCell_inc;
    else
        dEps = 1/dEps; % Change direction
    end 
    % - Run in that direcion    
    while(1)
%         [vCell_inc,dCell_inc] = evalSMH(x_sMH,ang_freqs,dEps*W_eps,nTheta,[],t_Thr,0);  
        [vCell_inc,dCell_inc] = evalSMH_sparse(x_sMH,ang_freqs,dEps*W_eps,nTheta,maxEigIdx,nn_cv,0);
        empLkh_inc = getCVlogL_gpu(x,ang_freqs,nTheta,trainSize,nv,vCell_inc,dCell_inc,lambda_c); 
        if empLkh_inc>=empLkh(cnt)
            first_W_eps_update = 0;
            cnt = cnt + 1;
            empLkh(cnt) = empLkh_inc;
            W_eps = dEps*W_eps;
            %disp(['Iteration ',num2str(cnt),': lambda_c = ',num2str(lambda_c),', epsilon = ',num2str(W_eps), ', Log-likelihood = ',num2str(empLkh(cnt))]);
            log_message('Iteration %d: lambda_c = %4.2f, epsilon = %4.2f, Log-likelihood = %4.2f',cnt,lambda_c,W_eps,empLkh(cnt));
            vCell = vCell_inc;
            dCell = dCell_inc;
        else
            if (first_W_eps_update==1)
                endFlag =1;
            end
            break;
        end
    end
    
    if (endFlag==1)
        break;
    end
end

%% Evaluate the steerable manifold harmonics
nn = ceil(nn_p*nImSMH/100);
if useWienerFilt
    x_sMH = PSWF_coeff_denSVS(:,1:nImSMH);
%    W_eps = sqrt(sum(wTot.^4))*nv;
else
    x_sMH = PSWF_coeff(:,1:nImSMH);
%     W_eps = sqrt(sum(rank))*nv;
end

% nThetaThr = 1e-2;

% - Find optimal parameters W_eps and nTheta
% [nTheta,~] = findOptEpsNTheta(x,ang_freqs,W_eps,nThetaThr);

% [vCell,dCell] = evalSMHOutOfMem_v2(x_sMH,ang_freqs,W_eps,nTheta,maxEigIdx,t_Thr);
% [vCell,dCell] = evalSMHOutOfMem_v3_sparse(x_sMH,ang_freqs,W_eps,nTheta,maxEigIdx,nn);
[vCell,dCell] = evalSMH_sparse(x_sMH,ang_freqs,W_eps,nTheta,maxEigIdx,nn,1);

%% Sort the steerable manifold harmonics according to their frequencies
dMat = [];
for i=0:max(ang_freqs)
    dMat = [dMat dCell{i+1}];
end
freqMat = repmat(0:max(ang_freqs),numel(dCell{1}),1);
[lambda,dSortIdx] = sort(dMat(:),'ascend');
freqMat = freqMat(:); freqMat = freqMat(dSortIdx);

evIdxMat = zeros(numel(lambda),max(ang_freqs)+1);
evIdxMat(1,freqMat(1)+1) = 1;
for i=2:numel(freqMat)
    evIdxMat(i,:) = evIdxMat(i-1,:);
    evIdxMat(i,freqMat(i)+1) = evIdxMat(i,freqMat(i)+1) + 1;
end

%% Extend to all data points by Nystrom's method
% vCell = smhExtNystrom(PSWF_coeff_denSVS,vCell,dCell,nImSMH,ang_freqs,W_eps,nTheta,normalizeDensity,denMat);
% vCellExt = smhExtNystrom(PSWF_coeff_denSVS,vCell,dCell,nImSMH,ang_freqs,W_eps,nTheta,normalizeDensity,0);
% vCell = vCellExt;

%% Find optimal K by a cross validation procedure
dl = 0.025;   % Resolution of lambda_c (cut-off frequnecy) in grid search between 0 and 1
% trainSize = nImages - 5e3;
trainSize = min(floor(nImages/2),nImSMH);
% trainSize = nImSMH;
% lambdaOpt = findCutoffCrossVal(PSWF_coeff,ang_freqs,vCell,dCell,dl,nTheta,nv,trainSize);
lambdaOpt = findCutoffCrossVal(PSWF_coeff(:,1:2*trainSize),ang_freqs,vCell,dCell,dl,nTheta,nv,trainSize);

%% Expanding the images by the steerable manifold harmonics, and performing sigular-value shrinkage on the expansion coefficients
% kOpt = size(evIdxMat,1);
% lambdaOpt = 0.9;
kOpt = find(lambda>lambdaOpt,1,'first')-1;
% [PSWF_coeff_denMan,~] = denoiseCoeffBySMH(PSWF_coeff,ang_freqs,evIdxMat(kOpt,:),vCell,nv,1:nImages);
log_message('Denoising by sMH regression...')
[PSWF_coeff_denMan,~] = denoiseCoeffBySMH(PSWF_coeff,ang_freqs,evIdxMat(kOpt,:),vCell,nv,1:nImSMH);

%% EM coeff estimation
log_message('Refining denoising by Expectation-Maximization...');
[ PSWF_coeff_denMan_EM,Likelihood,lScore] = denoiseCoeffBySMH_EM( PSWF_coeff,PSWF_coeff_denMan(:,1:nImSMH),ang_freqs,evIdxMat(kOpt,:),vCell,nv,2,256 );

%% View preliminary de-noising results for debugging
% testIdx = 5;
% I_denoise_man = mu + Psi(:,ang_freqs==0)*PSWF_coeff_denSVS(ang_freqs==0,testIdx) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff_denSVS(ang_freqs~=0,testIdx));
% image_d_disp = zeros(2*L,2*L);
% image_d_disp(r<=r_max) = real(I_denoise_man);
% figure; imagesc(image_d_disp); colormap(gray)
% 
% I_denoise_man = mu + Psi(:,ang_freqs==0)*PSWF_coeff_denMan(ang_freqs==0,testIdx) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff_denMan(ang_freqs~=0,testIdx));
% image_d_disp = zeros(2*L,2*L);
% image_d_disp(r<=r_max) = real(I_denoise_man);
% figure; imagesc(image_d_disp); colormap(gray)
% 
% I_denoise_man = mu + Psi(:,ang_freqs==0)*PSWF_coeff_denMan_EM(ang_freqs==0,testIdx) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff_denMan_EM(ang_freqs~=0,testIdx));
% image_d_disp = zeros(2*L,2*L);
% image_d_disp(r<=r_max) = real(I_denoise_man);
% figure; imagesc(image_d_disp); colormap(gray)

%% Form denoised images
I_denoise_svs = mu*ones(1,nImages) + Psi(:,ang_freqs==0)*PSWF_coeff_denSVS(ang_freqs==0,:) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff_denSVS(ang_freqs~=0,:));
projections_svs_den = zeros(2*L,2*L,nImages);
for i=1:nImages
    tmpIm = zeros(2*L,2*L);
    tmpIm(r<=r_max) = real(I_denoise_svs(:,i));
    projections_svs_den(:,:,i) = tmpIm;
end

I_denoise_man = mu*ones(1,size(PSWF_coeff_denMan,2)) + Psi(:,ang_freqs==0)*PSWF_coeff_denMan(ang_freqs==0,:) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff_denMan(ang_freqs~=0,:));
projections_man_den = zeros(2*L,2*L,size(I_denoise_man,2));
for i=1:size(I_denoise_man,2)
    tmpIm = zeros(2*L,2*L);
    tmpIm(r<=r_max) = real(I_denoise_man(:,i));
    projections_man_den(:,:,i) = tmpIm;
end

I_denoise_sPca = mu*ones(1,nImages) + Psi(:,ang_freqs==0)*PSWF_coeff(ang_freqs==0,:) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff(ang_freqs~=0,:));
projections_sPca = zeros(2*L,2*L,size(I_denoise_sPca,2));
for i=1:size(I_denoise_sPca,2)
    tmpIm = zeros(2*L,2*L);
    tmpIm(r<=r_max) = real(I_denoise_sPca(:,i));
    projections_sPca(:,:,i) = tmpIm;
end

I_denoise_man_EM = mu*ones(1,size(PSWF_coeff_denMan_EM,2)) + Psi(:,ang_freqs==0)*PSWF_coeff_denMan_EM(ang_freqs==0,:) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff_denMan_EM(ang_freqs~=0,:));
projections_man_EM_den = zeros(2*L,2*L,size(I_denoise_man_EM,2));
for i=1:size(I_denoise_man_EM,2)
    tmpIm = zeros(2*L,2*L);
    tmpIm(r<=r_max) = real(I_denoise_man_EM(:,i));
    projections_man_EM_den(:,:,i) = tmpIm;
end

%% Reshuffle back the images
% if shuffle
%     projections_man_EM_den(:,:,rndPermIdx) = projections_man_EM_den;
%     projections_man_den(:,:,rndPermIdx) = projections_man_den;
%     projections_svs_den(:,:,rndPermIdx) = projections_svs_den;
%     projections_sPca(:,:,rndPermIdx) = projections_sPca;
% end
