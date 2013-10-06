% Estimate the 1D isotropic autocorrelation function 
%
% Syntax: [R,x,cnt]=epsdR(vol,samples_idx,max_d,biasflag,verbose)
%
% Inputs:
%    vol - a tack of square pxp projections
%    samples_idx - is a list of pixel indices to use for autocorrelation estimation
%    biasflag - flag controls the bias of the estimate:
%    0    do not remove the mean from the samples before computing the
%         autocorrelation.
%    1    remove the mean
%    2    do not remove the mean but set the variance (autocorrelation at
%         distance 0) to be 1 (default).
%
%
% Outputs:
%    R - ACF samples
%    x - lag points
%    cnt - number of pixels available in each lag
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:  Yoel Shkolnisky, Mor Cohen
%________________________________________

% email: morcohe5@mail.tau.ac.il
% May 2011; Last revision: 04-May-2011
%
% See also: 

function [R,x]=epsdR_speedup(vol,samples_idx, biasflag, verbose)

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

if ~exist('biasflag','var')
    biasflag=2;
end

if ~exist('verbose','var')
    verbose=0;
end

% compute the number of elements used to compute each auto correlation
% entry. That is, in the sum R(k1,k2)=sum_{i,j} X_{i,j} X_{i+k1,j+k2},
% compute how many entries are in the sum for each (k1,k2).
% This is done by sertting the participating image samples to 1 and
% computing autocorrelation again.

if verbose
    fprintf('processing projections...');
end

vol=reshape(vol, p^2, K);
vol=vol(samples_idx, :);

if biasflag==1 % Remove mean
    vol = vol - repmat(mean(vol, 2), 1, K);
end

if biasflag==2 % Remove mean and normalize variance
    vol = vol - repmat(mean(vol, 2), 1, K);
    for i=1:K
        vol(:, i)=vol(:, i)/norm(vol(:, i));
    end;
end;

[I, J]=ind2sub([p, p], samples_idx);

C=1/K*(vol*vol');
l=length(I);
Dist = -2*(I*I') - 2*(J*J') + repmat(I.^2+J.^2, 1, l)+repmat((I.^2+J.^2).', l, 1); 
[dist]=unique(Dist(:));
R=zeros(length(dist), 1);
for i=1:length(dist)
    R(i)=sum(C(Dist(:)==dist(i)))/length(find(Dist(:)==dist(i)));
end;
x=sqrt(dist);

if verbose
    fprintf('done\n');
end


% % % %% reference code compute correlations the slow way and compare
% % %           
% % % % Count how many different correlation distances we have
% % % [I,J]=meshgrid(0:max_d,0:max_d);
% % % dsquare2=I.*I+J.*J;
% % % dsquare2=sort(unique(dsquare2(dsquare2<=max_d*max_d)));
% % % 
% % % corrs2=zeros(size(dsquare2,1),3);
% % % corrs2(:,1)=dsquare2;
% % %           
% % % ns=zeros(p,frame_size);
% % % n=numel(ns);
% % % samples=zeros(n,3);
% % % 
% % % tic
% % % for k=1:K
% % %     ns(:,:)=vol(:,1:frame_size,k);
% % %     [I,J]=ind2sub([p,frame_size],1:numel(ns));
% % %     samples(:,1)=I(:);
% % %     samples(:,2)=J(:);
% % %     samples(:,3)=ns(:);
% % %     samples(:,3)=samples(:,3)-mean(samples(:,3));
% % %     
% % %     % compute correlation of the samples
% % %     for j1=1:n        
% % %         for j2=j1:n
% % %             d2=(samples(j1,1)-samples(j2,1)).^2+...
% % %                 (samples(j1,2)-samples(j2,2)).^2;
% % % 
% % %             % take only positive shifts k1,k2 > 0. That is, the second
% % %             % sample in the product  X_{i,j} X_{i+k1,j+k2} should come from
% % %             % a pixel that is further down and right.
% % %             if (d2<=max_d*max_d) && (samples(j2,1)>=samples(j1,1)) && (samples(j2,2)>=samples(j1,2))
% % %                 v=samples(j1,3)*samples(j2,3);
% % %                 idx=find(corrs2(:,1)==d2);
% % %                 corrs2(idx,2)=corrs2(idx,2)+v;
% % %                 corrs2(idx,3)=corrs2(idx,3)+1;
% % %             end            
% % %         end
% % %         
% % %     end
% % %     fprintf('k=%d/%d\n',k, K);
% % % end
% % % 
% % % corrs2(corrs2(:,3)==0,:)=[];
% % % R2=corrs2(:,2)./corrs2(:,3);
% % % toc

