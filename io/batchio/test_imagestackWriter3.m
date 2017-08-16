% Large scale test of imagestackWriter
%
% Write a large number of images to the disk using imagestackWriter, and
% compute their hash on-the-fly. Once done. Read the images using
% imagestackReader, and compute the hash of the read images.
% The two hash codes should be the same.
%
% Yoel Shkolnisky, May 2016.

%% Write many images to disk
initstate;
n=65;
K=10000;
SNR=1; % Dummay SNR
blocksize=1000;
nimages=0;
hash1=0;
fname='tmp.mrc';
outstack=imagestackWriter(fname,K,1,102); % Cache size that does not divide K.


% Create a huge data set
Nblocks=round(K/blocksize);
for k=1:Nblocks
    % Generate the next block of images
    fprintf('Writing block of images %3d/%d\n',k,Nblocks);
    projs=cryo_gen_projections(n,blocksize,SNR);
    projs=single(projs);    
    nimages=nimages+blocksize;
    outstack.append(projs);
    
    % Update md5 of the stream
    opt.Input='bin';
    opt.Format='double';
    opt.Method='MD5';
    buf=[projs(:);hash1(:)];
    hash1=DataHash(buf,opt);

end

outstack.close;

fprintf('Written %d image to %s\n',nimages,fname);
fprintf('Hash code of written data: ');
hashhex1=dec2hex(hash1);
for k=1:16 
    fprintf('%s%s ',hashhex1(k,1),hashhex1(k,2));
end
fprintf('\n');


%% Read back the images
stack=imagestackReader(fname,100);
hash2=0;
projs=zeros(n,n,blocksize);

idx=1;
for k=1:Nblocks
    fprintf('Reading block of images %3d/%d\n',k,Nblocks);
    for j=1:blocksize
        projs(:,:,j)=stack.getImage(idx);
        idx=idx+1;
    end
    
    projs=single(projs);
    opt.Input='bin';
    opt.Format='double';
    opt.Method='MD5';
    buf=[projs(:);hash2(:)];
    hash2=DataHash(buf,opt);
    
end

fprintf('Read %d image from %s\n',idx-1,fname);
fprintf('Hash code of read data:    ');
hashhex2=dec2hex(hash2);
for k=1:16 
    fprintf('%s%s ',hashhex2(k,1),hashhex2(k,2));
end
fprintf('\n');

if hash1==hash2
    fprintf('test ok\n');
else
    fprintf('test FAILED\n');
end