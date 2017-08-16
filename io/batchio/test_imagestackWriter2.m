% Basic test of imagestackWriter
%
% Similar to test)imagestackWriter1, but write more images.
%
% Yoel Shkolnisky, May 2016.

% Generate projections.
initstate;
n=65;
K=500;
SNR=1; % Dummay SNR
projs=cryo_gen_projections(n,K,SNR);
projs=single(projs);

% Write to disk.
outstack=imagestackWriter('tmp.mrc',K,1,102); % Cache size that does not divide K.
outstack.append(projs);
outstack.close;

% Read saved images from MRC
projs2=ReadMRC('tmp.mrc');

% Compare. projs and projs2 should agree to the bit.
err=norm(projs(:)-projs2(:));
fprintf('err=%e\n',err);
if err==0
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end