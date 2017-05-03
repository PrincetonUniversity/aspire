function fname=gen_spca_simulation_data(Nprojs,n,chunksize)
%
% GEN_SPCA_SIMULATION_DATA  Generate clean data in MRC format.
%
% This function was originally developed to benchmark SPCA algorithms, but
% can be used whenever stacks of clean projections are needed.
%
% gen_spca_simulation_data(Nprojs,n,chunksize)
%   Generate an MRC file named projs_n.mrc containing Nprojs clean
%   projection of size nxn of the density map EMD2275. At any given moment,
%   not more than chunksize images are saved in memory.
%   Defaults: Nprojs=100, n=129, chunksize=10000
%
% Example:
%   gen_spca_simulation_data(10000,129)
%
% Yoel Shkolnisky, September 2016.

if ~exist('Nprojs','var')
    Nprojs=100; % Number of projections
end

if ~exist('n','var')
    n=129;  % Each projection is nxn
end

if ~exist('chunksize','var')
    chunksize=10000;  % Maximal numbers of images to keep in memory
end



% Read a density map.
% Ribosome structures to near-atomic resolution from thirty thousand
% cryo-EM particles
mapname=cryo_fetch_emdID(2275);
map=ReadMRC(mapname);
map=single(map);

log_message('Resampling map to %dx%dx%d',n,n,n);
map=cryo_downsample(map,[n n n]);


% Generate projections of size 257x257
% Generate projections at full resolution.
if mod(size(map,1),2)==0 % Make odd-sized
    map=map(1:end-1,1:end-1,1:end-1);
end

q=qrand(Nprojs);  % Generate Nprojs projections to orient.
log_message('Generating %d projections of size %dx%d',Nprojs,n,n);

%debugprojs=zeros(n,n,Nprojs);
fname=sprintf('projs_%d.mrc',n);
outstack=imagestackWriter(fname,Nprojs,1,chunksize); 

idx=0; % How many projectioned were generated so far.
while idx<Nprojs
    tmpk=min(chunksize,Nprojs-idx); % Generate images in chuncks of chunksize.
    log_message('Generating projections %d to %d',idx+1,idx+tmpk);
    projs=cryo_project(map,q(:,idx+1:idx+tmpk),n);
    projs=permute(projs,[2,1,3]);
    outstack.append(projs);
    %debugprojs(:,:,idx+1:idx+tmpk)=projs;
    idx=idx+tmpk;
end
outstack.close;
log_message('Saved %d images of size %dx%d to MRC file %s',Nprojs,n,n,fname);

% Compare images saved on disk to the ones in memory.
% proj=ReadMRC(fname);
% err=norm(proj(:)-debugprojs(:))/norm(debugprojs(:));
% log_message('Relative difference between in-memory and on-disk images %e',err);


% Downsample to 129x129