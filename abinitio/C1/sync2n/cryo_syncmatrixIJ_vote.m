function R22=cryo_syncmatrixIJ_vote(clmatrix,k1,k2,K3,L,rots_ref,is_perturbed)
%
% Compute the (k1,k2) rotation block of the cryo-EM syncronization matrix,
% using the list of images specified in K3. For example,
%   R11=cryo_syncmatrixIJ(clmatrix,1,1,1:K,360,rots_ref)
%
% rots_ref (optional) are the rotations used to computed the common lines
% matrix.
%
% Yoel Shkolnisky, August 2010.
%
% Revisions:
% Y.S. March 01, 2015  Optimize the code by replacing the loop that call
%       rotratio_eulerangle with a call to the vectorized function
%       rotratio_eulerangle_vec.

ref=0;
if exist('rots_ref','var')
    if rots_ref==0
        ref=0;
    else
        ref=1;  % Reference rotations are given.
    end
end


if ~exist('is_perturbed','var')
    is_perturbed=0;
end

% % Check the input
% sz=size(clmatrix);
% if numel(sz)~=2
%     error('clmatrix must be a square matrix');
% end
% if sz(1)~=sz(2)
%     error('clmatrix must be a square matrix');
% end
% 
% K=sz(1);
TOL=1.0E-12; % This tolerance is relevant only if L is very large (1.0e15),
% to check for accuracy. Otherwise, ignore it and the
% resulting output.

J=[1 0 0; 0 1 0; 0 0 -1]; % Reflection matrix

[goodK3,~,~]=cryo_voteIJ(clmatrix,L,k1,k2,K3,rots_ref,is_perturbed);


% %%% START reference code
% % The following code has bbenoptimized below
% Rks=zeros(9,numel(K3));
% kk=1;
% for count=1:numel(goodK3) % For some reason "for k3=goodK3" does not work!?
%     k3=goodK3(count);
%     if (k3~=k1) && (k3~=k2)
%         clk=clmatrix([k1 k2 k3],[k1 k2 k3]);
%         Rk=rotratio_eulerangle(clk,L);
%         
%         % rotratio_eulerangle returns a numeric error code in case of an error.
%         % Check that no error has occured.
%         if numel(Rk)==9
%             R=(Rk+J*Rk*J)/2;
%             Rks(:,kk)=R(:);
%             kk=kk+1;
%             
%             if ref
%                 % Compare resulting rotation computed using common
%                 % lines to the rotations computed using the true
%                 % rotations.
%                 R1ref=rots_ref(:,:,k1);
%                 R2ref=rots_ref(:,:,k2);
%                 inv_R1ref=R1ref.';
%                 inv_R2ref=R2ref.';
%                 Rref=(inv_R1ref.'*inv_R2ref+J*inv_R1ref.'*inv_R2ref*J)/2;
%                 err=norm(Rref-R,'fro')/norm(Rref,'fro');
%                 if (err>TOL) && (L>1.0e14)
%                     %%%fprintf('Wrong rotation: [k1=%d  k2=%d  k3=%d] err=%e  tol=%e\n',k1,k2,k3,err,TOL);
%                 end
%             end
%         end
%     end
% end
% 
% % Merge all rotations computed by all triplets (k1,k2,k3) into a
% % single rotation that takes the frame of image k1 into the frame
% % of image k2 .
% if kk>1
%     Rks=Rks(:,1:kk-1);
%     diff=Rks-repmat(mean(Rks,2),1,kk-1);
%     err=norm(diff(:))/norm(Rks(:));
%     if err>TOL
%         %%%                fprintf('Inconsistent rotations: [k1=%d  k2=%d] err=%e  tol=%e\n',k1,k2,err,TOL);
%     end
%     %            assert(norm(diff(:))/norm(Rks(:))<TOL);
%     Rk=reshape(mean(Rks,2),3,3);
%     % XXX Enforce rotation ???
% else
%     % Images k1 and k2 correspond to the same viewing direction and
%     % differ only by in-plane rotation. No triangle can be formed by
%     % any image k3. Thus, we find the (in-plane) rotation matrix
%     % between k1 and k2 by aligning the two images. The matrix is
%     % [ cos(theta) sin(theta) 0 ;...
%     %  -sin(theta) cos(theta) 0 ;...
%     %       0           0     1 ]
%     % Here we cheat and compute it using the rotations.
%     % If rots_ref is not given, just put zero.
%     Rk=zeros(3,3);
%     if ref
%         R1ref=rots_ref(:,:,k1);
%         R2ref=rots_ref(:,:,k2);
%         inv_R1ref=R1ref.';
%         inv_R2ref=R2ref.';
%         
%         % Multiply by 2 so that the upper-left 2x2 block of Rk is equal
%         % to the upper-left 2x2 block of
%         % inv_R1ref.'*inv_R2ref+J*inv_R1ref.'*inv_R2ref*J.
%         Rk=2.*inv_R1ref.'*inv_R2ref;
%     end
% end
% 
% Told=toc;
% %%%END reference code

% S is a 2Kx2K matrix, containing KxK blocks of size 2x2.
% The (i,j) block is given by [r11 r12; r12 r22], where
% r_{kl}=<R_{i}^{k},R_{j}^{l}>, k,l=1,2, namely, the dot product of
% column k of R_{i} and columns l of R_{j}. Thus, given the true
% rotations R_{1},...,R_{K}, S is decomposed as S=W^{T}W where
% W=(R_{1}^{1},R_{1}^{2},...,R_{K}^{1},R_{K}^{2}), where R_{j}^{k}
% is the k column of R_{j}.
% To extract the required block of S, we note that Rk is of the
% form
%  Rk= [ r11 r12 0;...
%        r21 r22 0;...
%         0   0 r33 ].
%

%%%% START optimization of the reference code above
%tic;

[Rs,goodratios]=rotratio_eulerangle_vec(clmatrix,k1,k2,goodK3,L);
Rk=sum(Rs,3)/size(Rs,3);

%Rks2=zeros(2,2,numel(goodratios));

if ref
    for kk=1:numel(goodratios) % For some reason "for k3=goodK3" does not work!?
        %    R=Rs(1:2,1:2,kk);
        %    Rks2(:,kk)=R(:);
        
        %    Rks2(:,:,kk)=Rs(1:2,1:2,kk);
        
        
        % Compare resulting rotation computed using common
        % lines to the rotations computed using the true
        % rotations.
        R1ref=rots_ref(:,:,k1);
        R2ref=rots_ref(:,:,k2);
        inv_R1ref=R1ref.';
        inv_R2ref=R2ref.';
        Rref=(inv_R1ref.'*inv_R2ref+J*inv_R1ref.'*inv_R2ref*J)/2;
        Rref=Rref(1:2,1:2);
        R=Rs(1:2,1:2,kk);
        err=norm(Rref-R,'fro')/norm(Rref,'fro');
        if (err>TOL) && (L>1.0e14)
            %%%fprintf('Wrong rotation: [k1=%d  k2=%d  k3=%d] err=%e  tol=%e\n',k1,k2,k3,err,TOL);
        end
    end
end

if numel(goodratios)>0
    %diff=Rks2-repmat(mean(Rks2,2),1,numel(goodratios));
    diff=Rs(1:2,1:2,:)-repmat(Rk(1:2,1:2),[1,1,numel(goodratios)]);
    tmpR=Rs(1:2,1:2,:);
    err=norm(diff(:))/norm(tmpR(:));
    if err>TOL
        %%%                fprintf('Inconsistent rotations: [k1=%d  k2=%d] err=%e  tol=%e\n',k1,k2,err,TOL);
    end
    %            assert(norm(diff(:))/norm(Rks(:))<TOL);
    %Rk=reshape(mean(Rks2,2),2,2);
    %Rk=mean(Rks2,3);
    % XXX Enforce rotation ???
else
    % Images k1 and k2 correspond to the same viewing direction and
    % differ only by in-plane rotation. No triangle can be formed by
    % any image k3. Thus, we find the (in-plane) rotation matrix
    % between k1 and k2 by aligning the two images. The matrix is
    % [ cos(theta) sin(theta) 0 ;...
    %  -sin(theta) cos(theta) 0 ;...
    %       0           0     1 ]
    % Here we cheat and compute it using the rotations.
    % If rots_ref is not given, just put zero.
    Rk=zeros(3,3);
    if ref
        R1ref=rots_ref(:,:,k1);
        R2ref=rots_ref(:,:,k2);
        inv_R1ref=R1ref.';
        inv_R2ref=R2ref.';
        
        % Multiply by 2 so that the upper-left 2x2 block of Rk is equal
        % to the upper-left 2x2 block of
        % inv_R1ref.'*inv_R2ref+J*inv_R1ref.'*inv_R2ref*J.
        Rk=2.*inv_R1ref.'*inv_R2ref;
    end
end
%Tnew=toc;

%%%%% END optimization

%fprintf('Told/Tnew=%5.3f\n',Told/Tnew);
%assert(norm(Rk-Rk2)==0);
R22=Rk(1:2,1:2);
