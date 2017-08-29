function S=cryo_syncmatrix(clmatrix,L,rots_ref)
%
% Construct the CryoEM synchronization matrix, given a common lines matrix
% clmatrix, that was constructed using angular resolution of L radial lines
% per image.
%
% rots_ref (optional) are the rotations used to computed the common lines
% matrix. 
%
% Yoel Shkolnisky, August 2010.

ref=0;
if exist('rots_ref','var')
    ref=1;  % Reference quaternions are given.
end

% Check the input 
sz=size(clmatrix);
if numel(sz)~=2
    error('clmatrix must be a square matrix');
end
if sz(1)~=sz(2)
    error('clmatrix must be a square matrix');
end

K=sz(1);
TOL=1.0E-12; % This tolerance is relevant only if L is very large (1.0e15), 
             % to check for accuracy. Otherwise, ignore it and the
             % resulting output. 
J=[1 0 0; 0 1 0; 0 0 -1]; % Reflection matrix

S=eye(2*K);

for k1=1:K-1
    for k2=k1+1:K
        
        kk=1; % Each triplet (k1,k2,k3) defines some rotation matrix from 
              % image k1 to image k2. kk counts the number of such
              % rotations, that is, the number of thrid images k3 that give
              % rise to a rotation from image k1 to image k2. Not all
              % images k3 give rise to such a rotation, since, for example
              % the resuling triangle can be too small.
              
        % We store all rotation matrices for testing purposes.
        % This is not really needed, as the matrices can be processed on
        % the fly. 
        Rks=zeros(9,K); 
                
        for k3=1:K                        
            if (k3~=k1) && (k3~=k2)
                clk=clmatrix([k1 k2 k3],[k1 k2 k3]);
                Rk=rotratio_eulerangle(clk,L);
                
                % Test that the above line is the same as the old rotratio:
%                 Rk2=rotratio(clk,L);
%                 err1=norm(Rk-Rk2,'fro');
%                 err2=norm(Rk-J*Rk2*J,'fro');
%                 assert(min(err1,err2)<1.0e-12,'err=%e',min(err1,err2));
                
                % rotratio_eulerangle returns a numeric error code in case of an error. 
                % Check that no error has occured.
                if numel(Rk)==9                    
                    R=(Rk+J*Rk*J)/2;
                    Rks(:,kk)=R(:);
                    kk=kk+1;
                    
                    if ref
                        % Compare resulting rotation computed using common
                        % lines to the rotations computed using the true
                        % rotations.
                        R1ref=rots_ref(:,:,k1);
                        R2ref=rots_ref(:,:,k2);
                        inv_R1ref=R1ref.';
                        inv_R2ref=R2ref.';   
                        Rref=(inv_R1ref.'*inv_R2ref+J*inv_R1ref.'*inv_R2ref*J)/2;
                        err=norm(Rref-R,'fro')/norm(Rref,'fro');
                        if (err>TOL) && (L>1.0e14)
                            fprintf('Wrong rotation: [k1=%d  k2=%d  k3=%d] err=%e  tol=%e\n',k1,k2,k3,err,TOL);
                        end
                    end
                end
            end
        end
        
        % Merge all rotations computed by all triplets (k1,k2,k3) into a
        % single rotation that takes the frame of image k1 into the frame
        % of image k2 .
        if kk>1
            Rks=Rks(:,1:kk-1);
            diff=Rks-repmat(mean(Rks,2),1,kk-1);
            err=norm(diff(:))/norm(Rks(:));
            if err>TOL
                fprintf('Inconsistent rotations: [k1=%d  k2=%d] err=%e  tol=%e\n',k1,k2,err,TOL);
            end
%            assert(norm(diff(:))/norm(Rks(:))<TOL);
            Rk=reshape(mean(Rks,2),3,3);
            % XXX Enforce rotation ???
        else
            % Images k1 and k2 correspond to the same viewing direction and
            % differ only by in-plane rotation. No triangle can be formed by
            % any image k3. Thus, we find the (in-plane) rotation matrix
            % between k1 and k2 by aligning the two images. The matrix is
            % [ cos(theta) sin(theta) 0 ;...
            %  -sin(theta) cos(theta) 0 ;...
            %       0           0     1 ]
            % Here we cheat and compute it using the quaternions.
            % If rots_ref is not given, just put zero.
            Rk=zeros(3,3);
            if ref
                R1ref=rots_ref(:,:,k1);
                R2ref=rots_ref(:,:,k2);
                inv_R1ref=R1ref.';
                inv_R2ref=R2ref.';
                
                % The upper-left 2x2 block of Rk is equal
                % to the upper-left 2x2 block of
                % (inv_R1ref.'*inv_R2ref+J*inv_R1ref.'*inv_R2ref*J)/2.
                Rk=inv_R1ref.'*inv_R2ref;
            end
        end
       
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
         R22=Rk(1:2,1:2);   
         
         % Put the block R22 in the right place in S.
         S(2*(k1-1)+1:2*k1,2*(k2-1)+1:2*k2)=R22;
         S(2*(k2-1)+1:2*k2,2*(k1-1)+1:2*k1)=R22.';
    end
end


