function [R,x,cnt]=cryo_epsdR(vol,samples_idx,max_d,verbose)
% CRYO_EPSDR  Estimate the 1D isotropic autocorrelation of an image stack.
%
% [R,x,cnt]=cryo_epsdR(vol,samples_idx,max_d,verbose)
%   Estimate the 1D isotropic autocorrelation function of a stack of
%   images. The samples to use in each image are given in samples_idx. The
%   correlation is computed up to a maximal distance of max_d.
%
% Inputs parameters:
%    vol          Stack of square pxp images.
%    samples_idx  List of pixel indices to use for autocorrelation
%                 estimation.
%    max_d        Correlations are computed up to a maximal distance of
%                 max_d pixels. Default p-1.
%    verbose   Set to nonzero to print progress messages. Default 0.
%
% Output parameters:
%    R    1D vector with samples of the isotropic autocorrelation function.
%    x    Distaces at which the samples of the autocorrelation function are
%         given. A vector of the same length as R.
%    cnt  Number of autocorrelation samples available for each distance.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Yoel Shkolnisky, October 2014.


p=size(vol,1);
if size(vol,2)~=p
    error('vol must be a stack of square images');
end

if ndims(vol)>3
    error('vol must be a 3D array');
end

K=1;
if ndims(vol)==3
    K=size(vol,3);
end

if ~exist('max_d','var')
    max_d=p-1; % Do not compute correlations between pixels that are more 
               % than max_d pixels apart
end
max_d=min(max_d,p-1);

if ~exist('verbose','var')
    verbose=0;
end

% Generate all possible squared distances. For a vertical shift i
% and horizontal shift j, dists(i,j) cotnains the corresponding
% isotropic correlation distance i^2+j^2. dsquare is then the set of all
% possible squared distances, that is, all distances that can be generated
% by integer steps in the horizontal and vertical directions.
[I,J]=meshgrid(0:max_d,0:max_d);
dists=I.*I+J.*J;
dsquare=sort(unique(dists(dists<=max_d*max_d)));

corrs=zeros(size(dsquare,1),1);   % corrs(i) containts the sum of all 
                                  % products of the form x(j)x(j+d), where
                                  % d=sqrt(dsquare(i)).
corrcount=zeros(size(dsquare,1),1); % Number of pairs x(j)x(j+d) for each d.
x=sqrt(dsquare); % Distances at which the correlations are computed.

% Create a distance maps whose value at the index (i,j) is the index in the
% array dsquare whose value is i^2+j^2. Pairs (i,j) that correspond to
% distance that are larger than max_d are indicated by (-1).
distmap=(-1)*ones(size(dists));
for i=0:max_d
    for j=0:max_d
        d=i*i+j*j;
        if d<=max_d*max_d
            idx=bsearch(dsquare,d-1.0e-13,d+1.0e-13);
            if isempty(idx) || numel(idx)>1
                error('Something went wrong');
            end
            distmap(i+1,j+1)=idx;
        end
    end
end
validdisits=find(distmap~=-1); % These are the distances (i,j) such that 
                               % i*i+j*j<=max_d*max_d.

% Compute the number of terms in the expression sum_{j}x(j)x(j+d) for each
% distance d. As the correlation is two-dimensioanl, we compute for each
% sum of the form  R(k1,k2)=sum_{i,j} X_{i,j} X_{i+k1,j+k2}, 
% how many summands are in the in it for each (k1,k2).
% This is done by setting the participating image samples to 1 and
% computing autocorrelation again.
mask=zeros(p);
mask(samples_idx)=1;
tmp=zeros(2*p+1);
tmp(1:p,1:p)=mask;
ftmp=fft2(tmp);
c=ifft2(ftmp.*conj(ftmp));
c=c(1:max_d+1,1:max_d+1);
c=round(c); % Now c contains the number of terms in the sum above for each 
            % two-dimensional distance (k1,k2). XXX rename c to Ncorr.


R=zeros(length(corrs),1); % Values of the isotropic autocorrelation 
                          % function. R(i) is the value of the ACF at
                          % distance x(i).

Ttot = 0; % Total elapsed time.

fprintf('Processing projections...');
Tlm = -10; % When the last messgae has been printed. Print progess messages no faster than every one second.
refreshmessage = 1; % Should the progess message be refreshed.

for k=1:K       
    
    tic;
    
    if verbose
        if (refreshmessage) || (k==K) % Always refresh on last image.
            if k==1
                processedimages=1;
            else
                processedimages=k-1;
            end
            
            str=sprintf('%5d/%5d\t ETA %2.0f (secs)',k,K,Ttot*(K-k)/processedimages);
            fprintf(str);
            Tlm=tic;
            refreshmessage=0;
        end
    end        
    
    proj=vol(:,:,k);

    % Mask all pixels that are not used to autocorrelation estimation.
    samples=zeros(p);
    samples(samples_idx)=proj(samples_idx);
        
    % Compute non-periodic autocorrelation of masked image with itself.
    tmp=zeros(2*p+1);
    tmp(1:p,1:p)=samples;
    ftmp=fft2(tmp);
    s=ifft2(ftmp.*conj(ftmp));
    s=s(1:max_d+1,1:max_d+1);  % Retain only relevant distances.
           
    % Accumulate all autocorrelation values R(k1,k2) such that
    % k1^2+k2^2=const (all autocorrelations of a certain distance).

    % Reference code:
    % This is the reference code. 
    % It was replaced ny the code below, which is about 6 times faster.
    % for j=1:length(dsquare)
    %     idx=find(dists==dsquare(j));
    %     corrs(j)=corrs(j)+sum(s(idx));
    %     corrcount(j)=corrcount(j)+sum(c(idx));
    % end
    
    % First optimization of the reference code.
    % It is about 4.5 times faster than the reference code.    
    % for i=1:max_d+1
    %    for j=1:max_d+1
    %         idx=distmap(i,j);
    %         if idx~=-1
    %             corrs(idx)=corrs(idx)+s(i,j);
    %             corrcount(idx)=corrcount(idx)+c(i,j);
    %         end
    %    end
    % end    
               
    for j=1:numel(validdisits)
        currdist=validdisits(j);
        dmidx=distmap(currdist);
        corrs(dmidx)=corrs(dmidx)+s(currdist);
        corrcount(dmidx)=corrcount(dmidx)+c(currdist);
    end
    
    Ttot=Ttot+toc;
    
    if verbose
        if (toc(Tlm)>1.0) || (k>=K-1) % Refresh every one second, or last image is the last one.
            refreshmessage=1;
            fprintf(repmat('\b',1,numel(str)));
        end
    end
end


% Remove zero correlation sums (distances for which we had no samples at
% that distance)
idx=find(corrcount~=0);
R(idx)=R(idx)+corrs(idx)./corrcount(idx);
cnt = corrcount(idx);

idx=find(corrcount==0);
R(idx)=[];
x(idx)=[];


if verbose
    fprintf('done\n');
end

