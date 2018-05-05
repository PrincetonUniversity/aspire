function [ projections_man_den,projections_svs_den,projections_sPca,projections_man_EM_den] = rotInvManDeNoise_v1( filePath,beta,T,N,nImSMH,maxEigIdx,useSPCA,useWienerFilt,shuffle,normalizeDensity,olN,nn_p)
%% Add path
addpath(genpath('/home/borisl/Documents/ManifoldDenoising/FPswfCoeffEval'));
addpath(genpath('/home/yoel/data/work/aspire/io/'));
addpath(genpath('/home/yoel/data/work/aspire/development/'));

%% Generate grids and sizes
stack=imagestackReader(filePath); 
sizeVec = stack.dim; 
sizeX = sizeVec(1); sizeY = sizeVec(2); sizeZ = sizeVec(3);

if ~isempty(N)
    nImages = min(N,sizeZ);
else
    nImages = sizeZ;
end
nImSMH = min(nImSMH,nImages);

L=floor(min(sizeX,sizeY)/2);
[x, y]=meshgrid(-L:L-1, -L:L-1);    % Even-numbered grid
% [x, y]=meshgrid(-L:L,-L:L); 

r=sqrt(x.^2+y.^2);
r_max = L;

%% Load data
sVarVec = zeros(1,nImages);
nVarVec = zeros(1,nImages);
projections = zeros(2*L,2*L,nImages);
for j = 1:nImages
        currImage=stack.getImage(j);
%         currImage = currImage(1:2*L,1:2*L)/sqrt(var(currImage(r>r_max)));
        projections(:,:,j) = currImage(1:2*L,1:2*L);
        sVarVec(j) = var(currImage(r<=r_max));
        nVarVec(j) = var(currImage(r>r_max));
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
end

%% Estimate noise variance and signal power
nv = mean(nVarVec);
sp = mean(sVarVec);

snr = 10*log10(sp/nv-1);

%% Map images to PSWF expansion coefficients
[ PSWF_coeff, Psi, ~, ang_freqs, ~, ~] = fastPswfCoeffEval( projections, x(r<=r_max)/r_max, y(r<=r_max)/r_max, r_max, beta, T, 64, 1 );

%% Perform Preliminary de-noising of PSWF coefficients by singular value shrinkage (for every angular index)
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
    clc; disp(['Performing preliminairy de-noising of angular index: ',num2str(m),', out of ',num2str(max(ang_freqs))]);
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
% I_denoise_pre = mu + Psi(:,ang_freqs==0)*PSWF_coeff_denSVS(ang_freqs==0,testIdx) + 2*real(Psi(:,ang_freqs~=0)*PSWF_coeff_denSVS(ang_freqs~=0,testIdx));
% image_d_disp = zeros(2*L,2*L);
% image_d_disp(r<=r_max) = real(I_denoise_pre);
% figure; imagesc(image_d_disp); colormap(gray)

%% Do mini-cross-validation on a subset of the data (over parameters lambda_c and epsilon)
nTheta = 256;
t_Thr = 1e-2;
% nMCV = nImSMH;
nMCV = 6e3;
nMCV = min(nMCV,nImages);
trainSize = round(nMCV/2);
nn_cv=ceil(nn_p*trainSize/100);

if useWienerFilt
    x_sMH = PSWF_coeff_denSVS(:,1:trainSize);
    x = PSWF_coeff(:,1:nMCV);
    W_eps = 4*sqrt(sum(wTot.^4))*nv;
else
    x_sMH = PSWF_coeff(:,1:trainSize);
    x = PSWF_coeff(:,1:nMCV);
    W_eps = 4*sqrt(sum(rank))*nv;
end

lambda_c = 0.5;

empLkh = [];
cnt = 1;
% [vCell,dCell,~] = evalSMH(x_sMH,ang_freqs,W_eps,nTheta,[],t_Thr,0);
[vCell,dCell] = evalSMHOutOfMem_v3_sparse(x_sMH,ang_freqs,W_eps,nTheta,maxEigIdx,nn_cv);
empLkh(cnt) = getCVlogL(x,ang_freqs,nTheta,trainSize,nv,vCell,dCell,lambda_c);

dlMCV = 0.05;
dEps = sqrt(2);
endFlag = 0;
while(1)        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Do lambda_c coordinate iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % - Determine direcion 
    empLkh_inc = getCVlogL(x,ang_freqs,nTheta,trainSize,nv,vCell,dCell,lambda_c+dlMCV); % Determine direcion
    if empLkh_inc>=empLkh(cnt)
        cnt = cnt + 1
        empLkh(cnt) = empLkh_inc;
        lambda_c = lambda_c + dlMCV;
    else
        dlMCV = -dlMCV; % Change direction
    end    
    % - Run in that direcion 
    while(1)
        empLkh_inc = getCVlogL(x,ang_freqs,nTheta,trainSize,nv,vCell,dCell,lambda_c+dlMCV);
        if empLkh_inc>=empLkh(cnt)
            cnt = cnt + 1
            empLkh(cnt) = empLkh_inc;
            lambda_c = lambda_c + dlMCV;
        else
            break;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Do W_eps coordinate iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    first_W_eps_update = 1;
    % - Determine direcion 
%     [vCell_inc,dCell_inc,~] = evalSMH(x_sMH,ang_freqs,dEps*W_eps,nTheta,[],t_Thr,0);  
    [vCell_inc,dCell_inc] = evalSMHOutOfMem_v3_sparse(x_sMH,ang_freqs,dEps*W_eps,nTheta,maxEigIdx,nn_cv);
    empLkh_inc = getCVlogL(x,ang_freqs,nTheta,trainSize,nv,vCell_inc,dCell_inc,lambda_c);     
    if empLkh_inc>=empLkh(cnt)
        first_W_eps_update = 0;
        cnt = cnt + 1
        empLkh(cnt) = empLkh_inc;
        W_eps = dEps*W_eps;
        vCell = vCell_inc;
        dCell = dCell_inc;
    else
        dEps = 1/dEps; % Change direction
    end 
    % - Run in that direcion    
    while(1)
%         [vCell_inc,dCell_inc] = evalSMH(x_sMH,ang_freqs,dEps*W_eps,nTheta,[],t_Thr,0);  
        [vCell_inc,dCell_inc] = evalSMHOutOfMem_v3_sparse(x_sMH,ang_freqs,dEps*W_eps,nTheta,maxEigIdx,nn_cv);
        empLkh_inc = getCVlogL(x,ang_freqs,nTheta,trainSize,nv,vCell_inc,dCell_inc,lambda_c); 
        if empLkh_inc>=empLkh(cnt)
            first_W_eps_update = 0;
            cnt = cnt + 1
            empLkh(cnt) = empLkh_inc;
            W_eps = dEps*W_eps;
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
%     W_eps = 4*sqrt(sum(wTot.^4))*nv;
else
    x_sMH = PSWF_coeff(:,1:nImSMH);
%     W_eps = 4*sqrt(sum(rank))*nv;
end

% nThetaThr = 1e-2;

% - Find optimal parameters W_eps and nTheta
% [nTheta,~] = findOptEpsNTheta(x,ang_freqs,W_eps,nThetaThr);

% [vCell,dCell] = evalSMHOutOfMem_v2(x_sMH,ang_freqs,W_eps,nTheta,maxEigIdx,t_Thr);
[vCell,dCell] = evalSMHOutOfMem_v3_sparse(x_sMH,ang_freqs,W_eps,nTheta,maxEigIdx,nn);

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
dl = 0.05;   % Resolution of lambda_c (cut-off frequnecy) in grid search between 0 and 1
% trainSize = nImages - 5e3;
trainSize = min(round(nImages/2),5e3);
lambdaOpt = findCutoffCrossVal(PSWF_coeff,ang_freqs,vCell,evIdxMat,lambda,dl,nTheta,nv,trainSize);

%% Expanding the images by the steerable manifold harmonics, and performing sigular-value shrinkage on the expansion coefficients
% kOpt = size(evIdxMat,1);
% lambdaOpt = 0.75;
kOpt = find(lambda>lambdaOpt,1,'first')-1;
% [PSWF_coeff_denMan,~] = denoiseCoeffBySMH(PSWF_coeff,ang_freqs,evIdxMat(kOpt,:),vCell,nv,1:nImages);
[PSWF_coeff_denMan,~] = denoiseCoeffBySMH(PSWF_coeff,ang_freqs,evIdxMat(kOpt,:),vCell,nv,1:nImSMH);

%% EM coeff estimation
% [ PSWF_coeff_denMan_EM,Likelihood] = denoiseCoeffBySMH_EM( PSWF_coeff,PSWF_coeff_denMan(:,1:nImSMH),ang_freqs,maxEigIdx,vCell_ext,nv,2,128 );
[ PSWF_coeff_denMan_EM,Likelihood] = denoiseCoeffBySMH_EM( PSWF_coeff,PSWF_coeff_denMan(:,1:nImSMH),ang_freqs,evIdxMat(kOpt,:),vCell,nv,2,128 );
% [ PSWF_coeff_denMan_EM,Likelihood] = denoiseCoeffBySMH_EM( PSWF_coeff,PSWF_coeff_denMan,ang_freqs,evIdxMat(kOpt,:),vCell,nv,2,128 );

%% View preliminary de-noising results for debugging
% testIdx = 1;
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
if shuffle
    projections_man_EM_den(:,:,rndPermIdx) = projections_man_EM_den;
    projections_man_den(:,:,rndPermIdx) = projections_man_den;
    projections_svs_den(:,:,rndPermIdx) = projections_svs_den;
    projections_sPca(:,:,rndPermIdx) = projections_sPca;
end