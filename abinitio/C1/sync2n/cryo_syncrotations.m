function [rotations,diff,mse,O]=cryo_syncrotations(S,rots_ref,verbose)
%
% Compute the rotations from the syncronization matrix S.
% rots_ref (optional) are the true rotations.
%
% Output parameters:
%   rotations   3x3xK array where rotations(:,:,k) is the k'th recovered
%               rotation.
%   diff        A vector of length K, whose k'th elements is the difference
%               between the recovered k'th rotations and the true k'th
%               rotation as given by rots_ref(:,:,k). The difference is
%               measured as the Frobenius norm between the k'th recovered
%               matrix and the k'th true matrix, after registering the
%               recovered rotations to the true ones.
%   mse         Mean squared error, computed as the mean of the squares of
%               the array diff.
%   verbose     Nonzero to print debug messages. Default is 0.
%
% Yoel Shkolnisky, August 2010.
%
% Revised Y.S December 2014     Add verbose flag.

TOL=1.0e-14;

if ~exist('verbose','var')
    verbose=0;
end

ref=0;
if exist('rots_ref','var')
    ref=1;  % Reference rotations are given.
end

% Check the input 
sz=size(S);
if numel(sz)~=2
    error('S must be a square matrix');
end
if sz(1)~=sz(2)
    error('S must be a square matrix');
end

if mod(sz(1),2)==1
    error('S must be a square matrix of size 2Kx2K');
end

K=sz(1)/2;

% S is supposed to be of rank 3. Extract the three eigenvectors
% corresponding to non-zero eigenvalues.
[V,D]=eigs(S,10);
[sorted_eigs, sort_idx] = sort(diag(D),'descend');

if verbose
    disp('Top eigenvalues: ');
    sorted_eigs(1:10)
end
V(:,1:3) = V(:,sort_idx(1:3));

% S is a 2Kx2K matrix, containing KxK blocks of size 2x2.
% The (i,j) block is given by [r11 r12; r12 r22], where
% r_{kl}=<R_{i}^{k},R_{j}^{l}>, k,l=1,2, namely, the dot product of
% column k of R_{i} and columns l of R_{j}. Thus, given the true
% rotations R_{1},...,R_{K}, S is decomposed as S=W^{T}W where
% W=(R_{1}^{1},R_{1}^{2},...,R_{K}^{1},R_{K}^{2}), where R_{j}^{k}
% is the k column of R_{j}. Therefore, S is a rank-3 matrix, and thus, it
% three eigenvectors that correspond to non-zero eigenvealues, are linear
% combinations of the column space of S, namely, W^{T}.

% According to the structure of W^{T} above, the odd rows of V, denoted V1,
% are a linear combination of the vectors R_{i}^{1}, i=1,...,K, that is of
% column 1 of all rotation matrices. Similarly, the even rows of V,
% denoted, V2, are linear combinations of R_{i}^{1}, i=1,...,K.
V1 = V(1:2:2*K,1:3);
V2 = V(2:2:2*K,1:3);

V1 = V1';
V2 = V2';

% We look for a linear transformation (3 x 3 matrix) A such that
% A*V1'=R1 and A*V2=R2 are the columns of the rotations matrices.
%
% Therefore:
%
% V1 * A'*A V1' = 1
% V2 * A'*A V2' = 1
% V1 * A'*A V2' = 0
%
% These are 3*K linear equations for to 9 matrix entries of A'*A
% Actually, there are only 6 unknown variables, because A'*A is symmetric.

% 3*K equations in 9 variables (3 x 3 matrix entries).
equations = zeros(3*K,9);

for i=1:3;
    for j=1:3;
        equations(1:3:3*K, 3*(i-1)+j) = V1(i,:) .* V1(j,:);
        equations(2:3:3*K, 3*(i-1)+j) = V2(i,:) .* V2(j,:);
        equations(3:3:3*K, 3*(i-1)+j) = V1(i,:) .* V2(j,:);
    end;
end;

% Truncate from 9 variables to 6 variables corresponding
% to the upper half of the matrix A'*A
truncated_equations = equations(:, [1, 2, 3, 5, 6, 9]);

% b = [1 1 0 1 1 0 ...]' is the right hand side vector
b = ones(3*K,1);
b(3:3:3*K) = zeros(K,1);

% Find the least squares approximation
ATA_vec = truncated_equations \ b;

% From the vectroized matrix ATA_vec construct back the matrix A'*A.
ATA = zeros(3);
ATA(1,1) = ATA_vec(1);
ATA(1,2) = ATA_vec(2);
ATA(1,3) = ATA_vec(3);
ATA(2,1) = ATA_vec(2);
ATA(2,2) = ATA_vec(4);
ATA(2,3) = ATA_vec(5);
ATA(3,1) = ATA_vec(3);
ATA(3,2) = ATA_vec(5);
ATA(3,3) = ATA_vec(6);

% The Cholesky decomposition of A'*A gives A
A = chol(ATA);

% Recover the rotations. The first two columns of all rotation matrices are
% given by unmixing V1 and V2 using A. The third column is the cross
% product of the first two.
rotations = zeros(3,3,K);
R1 = A*V1;
R2 = A*V2;
R3 = cross(R1,R2);

for k=1:K
    rotations(:,1,k)=R1(:,k);
    rotations(:,2,k)=R2(:,k);
    rotations(:,3,k)=R3(:,k);
end

% Make sure that we got rotations.
for k=1:K
    R=rotations(:,:,k);
    erro=norm(R*R.'-eye(3));
    if erro>TOL
%%%        fprintf('Trnaformation %d is not orthogonal, err=%e  tol=%e\n',k,erro,TOL);
    end
    
    errd=abs(det(R)-1);
    if errd>TOL
%%%        fprintf('Determinant of %d diffrs from 1, err=%e  tol=%e\n',k,errd,TOL);    
    end
    
    % Enforce R to be a rotation (in case the error is large)
    [U,~,V]=svd(R);
    rotations(:,:,k)=U*V.';
end


% If reference rotations are given, compare the resluting rotations against
% true ones.

if ref   
    J=[1 0 0; 0 1 0; 0 0 -1]; % Reflection matrix
    
    for k1=1:K-1
        for k2=k1+1:K
            R1=rotations(:,:,k1);
            R2=rotations(:,:,k2);
            R=R1.'*R2;
            
            R1ref=rots_ref(:,:,k1);
            R2ref=rots_ref(:,:,k2);
            inv_R1ref=R1ref.';
            inv_R2ref=R2ref.';
            
            % The resulting rotations should satisfy the same ratio
            % equaions as the true orientations. Specifically, the ration
            % of each pair of estimated rotation should equal one of the
            % following two rotations:
            Rref1=inv_R1ref.'*inv_R2ref;
            Rref2=J*Rref1*J;
            
            err1=norm(Rref1-R,'fro')/norm(Rref1,'fro');
            err2=norm(Rref2-R,'fro')/norm(Rref2,'fro');
            if (err1>TOL) && (err2>TOL)
%                fprintf('Ratio error [k1=%d  k2=%d] err=%e  tol=%e\n',...
%                    k1,k2,min(err1,err2),TOL);
            end
        end
    end
    
    % Register estimated rotations to true ones, and compute the difference
    % between the two.
    rot=zeros(3*K,3);  % The K estimated rotation matrices stacked as a matrix with 3K rows.
    rot1=zeros(3*K,3); % True true K rotation matrices stacked as a matrix with 3K rows.
    rot2=zeros(3*K,3); % Reflected matrices of rot1, which are also a valid solution for the rotations.
    
    for k=1:K
        R=rotations(:,:,k);
        rot(3*(k-1)+1:3*k,:)=R.';
        Rref=rots_ref(:,:,k);
        rot1(3*(k-1)+1:3*k,:)=Rref;
        rot2(3*(k-1)+1:3*k,:)=J*Rref*J;
    end
    
    % Compute the two possible orthogonal matrices which register the
    % estimated rotations to the true ones.
    O1=rot.'*rot1./K;       
    O2=rot.'*rot2./K;
    
    % We are registering one set of rotations (the estimated ones) to
    % another set of rotations (the true ones). Thus, the transformation
    % matrix between the two sets of rotations should be orthogonal. This
    % matrix is either O1 if we recover the non-reflected solution, or O2,
    % if we got the reflected one. In any case, one of them should be
    % orthogonal.
    
    if verbose
        err1=norm(O1*O1.'-eye(3));
        err2=norm(O2*O2.'-eye(3));
        if (err1>TOL) && (err2>TOL)
            fprintf('Registering matrix is not orthogonal, err=%e  tol=%e\n',...
                min(err1,err2),TOL);
        end
        
        
        errd1=abs(det(O1)-1);
        errd2=abs(det(O2)-1);
        if (errd1>TOL) && (errd2>TOL)
            fprintf('Determinant of registering matrix is not 1, err=%e  tol=%e\n',...
                min(errd1,errd2),TOL);
        end
    end
    
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
    for k=1:K
        R=rotations(:,:,k);
        Rref=rots_ref(:,:,k);
        if flag==2
            Rref=J*Rref*J;
        end
        diff(k)=norm(R.'*O-Rref,'fro');
        mse=mse+diff(k).^2;
    end
    mse=mse/K;
%    hist(diff)
    
end
