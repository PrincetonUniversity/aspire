function test_cryo_prewhiten_outofcore
%
% Test the function cryo_prewhiten_outofcore by comapring its output to
% cryo_prewhite. Both functions shouls have the same output.
%
% Yoel Shkolnisky, May 2016.

N=65; % Each noise image is of size NxN.
Ks=1000; % Number of noise images to generate.

initstate;

% Generate a stack of noise images
noise=noise_exp2d(N,max(Ks),1);
instackname=tempname;
WriteMRC(single(noise),1,instackname); % Save the projections to an MRC file.
noise=ReadMRC(instackname);
%noise=double(noise); % Although we read single precision numbers, cast to 
    % double precision since cryo_normalize_background_outofcore uses
    % double precision internally, and we want both
    % cryo_normalize_background and cryo_normalize_background_outofcore to
    % have exactly the same roundoff error.

% Estimate power spectrum of the noise.
psd = cryo_noise_estimation(noise);

% Prewhiten images in-memory
pw1=cryo_prewhiten(noise,psd);

% Prewhiten images out-of-memory
outstackname=tempname;
cryo_prewhiten_outofcore(instackname,outstackname,psd);
pw2=ReadMRC(outstackname);

err=norm(pw1(:)-pw2(:));  % Should be zero.
fprintf('diff between two methods = %5.3e\n',err);
if err<1.0e-10
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end
