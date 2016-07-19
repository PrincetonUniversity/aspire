function [regrot,mse,diff,O,flag]=register_rotations(rotations,rotations_ref)
%REGISTER_ROTATIONS Register estimated orientations to reference ones.
% [regrot,mse,diff,O]=REGISTER_ROTATIONS(rotations,rotations_ref) Finds the
% orthogonal transformation that best aligns the estimated rotations to the
% reference rotations. "regrot" are the estimated rotations after registering
% them to the reference ones, O is the optimal orthogonal 3x3 matrix to
% align the two sets, mse is given by
%      mse=sum_{i} \norm(O.'*rotations(i)-rotations_ref(i),'fro')^2
% and diff is the square root of each term.
%
% If flag==2 then J conjugacy is required.
% Yoel Shkolnisky, October 2013.

K=size(rotations,3);

if any(size(rotations)-size(rotations_ref))
    error('rotations and rotations_ref must have same dimensions');
end

J=[1 0 0; 0 1 0; 0 0 -1]; % Reflection matrix

O1=zeros(3,3);
O2=zeros(3,3);

for k=1:K
   R=rotations(:,:,k);
   Rref=rotations_ref(:,:,k);
   O1=O1+R*Rref.';
   O2=O2+(J*R*J)*Rref.';
end


% Compute the two possible orthogonal matrices which register the
% estimated rotations to the true ones.
O1=O1./K;
O2=O2./K;

% We are registering one set of rotations (the estimated ones) to
% another set of rotations (the true ones). Thus, the transformation
% matrix between the two sets of rotations should be orthogonal. This
% matrix is either O1 if we recover the non-reflected solution, or O2,
% if we got the reflected one. In any case, one of them should be
% orthogonal.

err1=norm(O1*O1.'-eye(3));
err2=norm(O2*O2.'-eye(3));

% In cany case, enforce the registering matrix O to be a rotation.
if err1<err2
    [U,~,V]=svd(O1); % Use O1 as the registering matrix
    flag=1;
else
    [U,~,V]=svd(O2); % Use O2 as the registering matrix
    flag=2;
end
O=U*V.';

% Plot estimation errors
diff=zeros(K,1);
mse=0;
regrot=zeros(size(rotations));
for k=1:K
    R=rotations(:,:,k);
    Rref=rotations_ref(:,:,k);
    if flag==2
        R=J*R*J;
    end
    regrot(:,:,k)=O.'*R;
    diff(k)=norm(O.'*R-Rref,'fro');
    mse=mse+diff(k).^2;
end
mse=mse/K;

 