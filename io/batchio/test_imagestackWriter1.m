% Basic test of imagestackWriter
%
% Generate a small number of images,save them using imagestackWriter, load
% them using ReadMRC and compare the original (in-memory) stack with the
% one that has been written-read.
% The difference should be zero (both arrays should agree to the bit).

% Generate projections.
initstate;
n=65;
K=5;
SNR=1; % Dummay SNR
projs=cryo_gen_projections(n,K,SNR);
projs=single(projs);

% Write to disk.
outstack=imagestackWriter('tmp.mrc',K,1,50);
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