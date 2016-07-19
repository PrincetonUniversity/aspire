function [Mkj,Ckj,Cjk]=commonline_R_vec(Rks,Rj,L,delta)
%
% COMMONLINE_R_VEC  Compute common lines from rotation matrices.
%
% [Mkj,Ckj,Cjk]=commonline_R_vec(Rks,Rj,L,delta)
%   Compute the common lines induced by the rotation matrix Rj with all
%   rotations matrices in the array Rks. The size of all returned arrays is
%   size(Rks,3)x1. Mkj is a mask array with 1 for each element in Rks whose
%   viewing direction has angle at least acos(delta) with the viewing
%   direction Rj (that is, the cosine between the
%   viewing directions is smaller than delta). Ckj is the common line
%   between each Rks and Rj in the coordinates of Rks. Cjk is the common
%   line between each Rks and Rj in the coordinates of Rj. For entries for
%   which Mkj==0, Ckj and Cjk contain -1.
%
% The returned common lines are 1-based and not zero based as
% commonlines_R.
%
% Yoel Shkolnisky, June 2016.

Nrots=size(Rks,3);

% Rks_vec=zeros(3,3*Nrots);
% Rks_t_vec=zeros(3*Nrots,3);
% for k=1:Nrots
%     Rks_vec(:,3*(k-1)+1:3*k)=Rks(:,:,k);
%     Rks_t_vec(3*(k-1)+1:3*k,:)=(Rks(:,:,k)).';
% end

% Rks_vec_1=reshape(Rks,3,3*Nrots);
% Rks_t_vec_1=(Rks_vec_1).';
% 
% assert(norm(Rks_vec-Rks_vec_1)==0);
% assert(norm(Rks_t_vec-Rks_t_vec_1)==0);

% The following two lines are an optimized version of the above block. Note
% that the above block also include the required test code.
Rks_vec=reshape(Rks,3,3*Nrots);
Rks_t_vec=(Rks_vec).';

Mkj=zeros(Nrots,1);     % Pairs of rotations that are not "too close"

Rk3s=Rks_vec(:,3:3:end);
% Here essentially starts a vectorized implementation of the
% function commonline_R.
Rj3=Rj(:,3);
clvec = [Rk3s(2,:).*Rj3(3)-Rk3s(3,:).*Rj3(2);
    Rk3s(3,:).*Rj3(1)-Rk3s(1,:).*Rj3(3);
    Rk3s(1,:).*Rj3(2)-Rk3s(2,:).*Rj3(1)];

% No need to normalize clvec as the normalization does not affect the
% atan2 below.

clvec2=zeros(3*Nrots,3);

clvec2(1:3:end,:)=clvec.';
clvec2(2:3:end,:)=clvec.';
clvec2(3:3:end,:)=clvec.';

ckjs=sum(bsxfun(@times,Rks_t_vec,clvec2),2);
assert(norm(ckjs(3:3:end))<1.0e-13)

cjks=Rj.'*clvec;
cjks=cjks.';
assert(norm(cjks(:,3))<1.0e-13)

alphakj=atan2(ckjs(2:3:end),ckjs(1:3:end));
alphajk=atan2(cjks(:,2),cjks(:,1));

PI=4*atan(1.0);
alphakj=alphakj+PI; % Shift from [-pi,pi] to [0,2*pi].
alphajk=alphajk+PI;

ckj=alphakj/(2*PI)*L;
cjk=alphajk/(2*PI)*L;

Ckj=mod(round(ckj),L)+1;
Cjk=mod(round(cjk),L)+1;

Rtk3=Rks_t_vec(:,3);
Rtk3=reshape(Rtk3,3,Nrots);
Rtj3=(Rj(3,:)).';
idx1=find(sum(bsxfun(@times,Rtk3,Rtj3))<delta);
Mkj(idx1)=1;

idx2=find(sum(bsxfun(@times,Rtk3,Rtj3))>=delta);
Ckj(idx2)=-1;
Cjk(idx2)=-1;
