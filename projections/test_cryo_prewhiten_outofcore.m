function test_cryo_prewhiten_outofcore
%
% Test the function cryo_prewhiten_outofcore by comapring its output to
% cryo_prewhite. Both functions shouls have the same output.
%
% Yoel Shkolnisky, May 2016.

N=128; % Each noise image is of size NxN.
Ks=1000; % Number of noise images to generate.

initstate;

% Generate a stack of noise images
noise=noise_exp2d(N,max(Ks),1);
instackname=tempname;
WriteMRC(single(noise),1,instackname); % Save the projections to an MRC file.
noise=ReadMRC(instackname);

% Estimate power spectrum of the noise.
psd = cryo_noise_estimation(noise);

% Prewhiten images in-memory
pw1=cryo_prewhiten(noise,psd);

% Prewhiten images out-of-memory
outstackname=tempname;
cryo_prewhiten_outofcore(instackname,outstackname,psd);
pw2=ReadMRC(outstackname);

% Check error
err=zeros(Ks,1);
for k=1:Ks
    err(k)=norm(pw1(:,:,k)-pw2(:,:,k))/norm(pw1(:,:,k));
end
merr=max(abs(err));
fprintf('max relative diff between two methods = %5.3e\n',merr);
if merr<5*eps(class(pw2))
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end
