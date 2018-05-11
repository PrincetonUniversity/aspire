function S=cryo_syncmatrix_vote_old(clmatrix,L,rots_ref,is_perturbed)
%
% Construct the CryoEM synchronization matrix, given a common lines matrix
% clmatrix, that was constructed using angular resolution of L radial lines
% per image.
%
% rots_ref (optional) are the rotations used to computed the common lines
% matrix. 
%
% is_perturbed descrbies which common lines are correct. It is of the same
% size as clmatrix, and the (i,j) entry is 0 is the common line between
% images i and j is correct and 1 otherwise.
%
% Like cryo_syncmatrix_p3, but uses volting to compute the IJ element.
% 
% Yoel Shkolnisky, September 2010.

if ~exist('rots_ref','var')
    rots_ref=0;
end

if ~exist('is_perturbed','var')
    is_perturbed=0;
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
S=eye(2*K);   

for k1=1:K-1
    %fprintf('Process image %d out of %d\n',k1,K-1);
    Stmp=zeros(2,2,K);    
    parfor k2=k1+1:K
    %for k2=k1+1:K
        Stmp(:,:,k2)=cryo_syncmatrixIJ_vote_old(clmatrix,k1,k2,1:K,L,rots_ref,is_perturbed);
    end
    
    for k2=k1+1:K
        R22=Stmp(:,:,k2);
        S(2*(k1-1)+1:2*k1,2*(k2-1)+1:2*k2)=R22;
        S(2*(k2-1)+1:2*k2,2*(k1-1)+1:2*k1)=R22.';
    end

end


