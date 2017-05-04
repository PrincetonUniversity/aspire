% Test the function cryo_pft_outofcore by comparing it to cryo_pft
%
% Yoel Shkolnisky, April 2017.

K=100;
n=129;

% Generate polar Fourier rays cryo_pft 
images=rand(n,n,K,'single');
n_theta=360;
n_r=ceil(n/2);
pf=cryo_pft(images,n_r,n_theta,'single');
pf=single(pf); % cryo_pft_outofcore returns single precision.

% Write Fourier rays into an MRC file
pfin='temppf.mrc';
imstack=imagestackWriterComplex(pfin,K);
imstack.append(pf);
imstack.close;

% Normalize Fourier rays in-core
pf_norm1=cryo_raynormalize(pf);

% Normalize Fourier rays out-of-core
pfout='temppfnorm.mrc';
cryo_raynormalize_outofcore(pfin,pfout);

% Compare outputs of in-core and out-of-core functions
pfstack=imagestackReaderComplex(pfout);
n_projs=pfstack.dim(3);
pf_norm2=pfstack.getImage(1:n_projs);

err=norm(pf_norm1(:)-pf_norm2(:));
fprintf('Values err=%e\n',err);
if err==0
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end
