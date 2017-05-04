% Test the function cryo_pft_outofcore by comparing it to cryo_pft
%
% Yoel Shkolnisky, April 2017.

K=100;
n=129;
images=rand(n,n,K,'single');

n_theta=360;
n_r=ceil(n/2);

% Compute cryo_pft directly
[pf1,freqs1]=cryo_pft(images,n_r,n_theta,'single');
pf1=single(pf1); % cryo_pft_outofcore returns single precision.

% Write images into an MRC file
imname='tempin.mrc';
imstack=imagestackWriter(imname,K);
imstack.append(images);
imstack.close;

% Compute out-of-core pft
pfname='temppf.mrc';
freqs2=cryo_pft_outofcore(imname,pfname,n_r,n_theta);

% Make sure both functions return the same frequnecies
err=norm(freqs1(:)-freqs2(:));
fprintf('freqs err=%e\n',err);
if err==0
    fprintf('freqs ok\n');
else
    fprintf('freqs FAILED\n');
end

% Compare outputs of polar Fourier transforms
pfstack=imagestackReaderComplex(pfname);
n_projs=pfstack.dim(3);
pf2=pfstack.getImage(1:n_projs);

err=norm(pf1(:)-pf2(:));
fprintf('Values err=%e\n',err);
if err==0
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end
