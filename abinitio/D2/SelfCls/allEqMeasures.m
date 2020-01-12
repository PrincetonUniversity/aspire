
%% Compute a measure of how much an image from data is close to be an equator
%Input:  L = Sampling resolution        
%        C = L by L/2 fourier ray correlations matrix
function [corrs_means,normals_corrs_max]=allEqMeasures(C,L)

% First compute the eq measure (corrs(scl-k,scl+k) for k=1:90)
% An eqautor image of a D2 molecule has the following property: If t_i is
% the angle of one of the rays of the self common line the all the pairs of
% rays of the form (t_i-k,t_i+k) for k=1:90 are identical. For each t_i we
% average over correlations between the lines (t_i-k,t_i+k) for k=1:90
% to measure the likelihood that the image is an equator and the ray (line)
% with angle t_i is a self common line. 
% (This first loop can be done once outside this function and then pass
% idx as an argument). 
c1=1;
idx=zeros(2,90,180);
idx_1=mod([(c1-(1:90));(c1+(1:90))],360);
idx(:,:,1)=idx_1;
for k=1:179
    idx(:,:,1+k)=idx_1+k;
end
idx=mod(idx,360);

%Take care of zero entries (=360)
idx1=squeeze(idx(1,:,:));
idx2=squeeze(idx(2,:,:));
idx1=idx1(:);
idx2=idx2(:);
zero_idx1=abs(idx1)<1e-7;
zero_idx2=abs(idx2)<1e-7;
idx1(zero_idx1)=360;
idx2(zero_idx2)=360;

%Make all Rj coordinates <=180 and compute linear indices for corrrelations
bigger_than_180=idx2>180;
idx1(bigger_than_180)=mod(idx1(bigger_than_180)+180,360);
zero_idx1=abs(idx1)<1e-7;
idx1(zero_idx1,1)=360;
idx2(bigger_than_180)=idx2(bigger_than_180)-180;
corrs_idx=sub2ind([L,L/2],idx1,idx2);

%Compute correlations
corrs=C(corrs_idx);
corrs=reshape(corrs,[90,180]);
tmp_corrs_means=mean(corrs,1);

% Now compute correlations for normals to scls. 
% An eqautor image of a D2 molecule has the additional following property:
% The normal line to a self common line in 2D Fourier plane is real values
% and both of its rays have identical values. We use the correlation
% between one Fourier ray of the normal to a self common line candidate t_i
% with its anti-podal as an additional way to measure if the image is an
% eqautor and t_i+0.5*pi is the normal to its self common line. 
r=2;
c3=c1+180;
% if c1<91+r
%     c3=c1+180;
% end
normal2scl1_idx=mod((c3-(90-r:90+r))',360); %(c3-(90-r:90+r))'+180]
normals2scls_idx=zeros(2*r+1,180);
normals2scls_idx(:,1)=normal2scl1_idx;
for k=1:179
   normals2scls_idx(:,1+k)=normal2scl1_idx+k;
end

%Take care of zero entries for normals to scls(=360)
normals2scls_idx=normals2scls_idx(:);
zero_idx=abs(normals2scls_idx)<1e-7;
normals2scls_idx(zero_idx)=360;
normals2scls_idx(normals2scls_idx<=180)=...
    normals2scls_idx(normals2scls_idx<=180)+180;
corrs_idx=sub2ind([L,L/2],normals2scls_idx,normals2scls_idx-180);
normals_corrs=reshape(C(corrs_idx),2*r+1,180);
normals_corrs_max=gather(max(normals_corrs,[],1));

%Compute self common lines correlations
% idx=[181:360;1:180];
% scls_lin_idx=sub2ind([360,180],idx(1,:),idx(2,:));
% scls_corrs=gather(C(scls_lin_idx));

%Finally compute top view corrs
%tv_corrs=repmat(0.5*(tmp_corrs_means(1:90)+tmp_corrs_means(91:180)),1,4);
%tv_corrs=gather(tv_corrs);

corrs_means=gather(tmp_corrs_means);


