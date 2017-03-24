function clerr=cryo_syncconsistency(rotations,clmatrix,L)
% Compute common lines consistency for given rotations and a given common
% lines matrix.
%
% The function takes rotations Ri, i=1,...,K, and computes for each pair of
% common lines cij and cji the angle (in degrees) between Ricij and Rjcji.
%
% Input parameters:
% rotations     Array of 3x3xK where rotations(:,:,i) is the rotation
%               matrix Ri.
% clmatrix      Common lines matrix, as computed, for example, by
%               cryo_clmatrix.
% L             Angular resolution - number of Fourier computed for each
%               image when searching for common lines.
%
% Output parameters:
% clerr         Consistency statistics. An array of K*(K-1)/2 row and three
%               columns. Each pair of images (i,j) has a row in the
%               table. The first column in the table contains i, the
%               second contains j, and the third contains the angles in
%               degrees between Ricij and Rjcji, where cij and cji are the
%               common line between images i and j.
%
% Yoel Shkolnisky, April 2012.
%
% Revision:
% Y.S. March 2015   Optimization for speed (about x10 faster now).


K=size(clmatrix,1);
clerr=zeros(K*(K-1)/2,3); % Common lines embedding errors for each of the
% K*(K-1)/2 common lines pairs. Each row of this array contains
% (i,j,e), where i and j are the indices of the common lines, and e is
% the dot product between R_{i}c_{ij} and R_{j}c_{ji}.
idx=0;

Rjs=sparse(3*K,3*K);
for kk=1:K
    Rjs(3*(kk-1)+1:3*kk,3*(kk-1)+1:3*kk)=rotations(:,:,kk);
end

for k1=1:K-1
    %fprintf('Processing k1=%d\n',k1);
    K2=k1+1:K;
    nnzidx=find(clmatrix(k1,K2)~=0); % Process only common lines that have
    % not been rejected in previous iterations (if there were any such
    % iterations), or during voting.
    alphai=2*pi*(clmatrix(k1,K2)-1)/L;
    alphai=alphai(:);
    alphaj=2*pi*(clmatrix(K2,k1)-1)/L;
    
    cij=zeros(numel(nnzidx),3);
    cij(:,1)=cos(alphai);
    cij(:,2)=sin(alphai);
    cij=cij.';
    
    cji=zeros(numel(nnzidx),3);
    cji(:,1)=cos(alphaj);
    cji(:,2)=sin(alphaj);
    cji=cji.';
    
    cjivec=zeros(3*K,1);
    cjivec(3*k1+1:end)=cji(:);
    cljvec=Rjs*cjivec;
    cljvec=cljvec(3*k1+1:end);
    cljvec=reshape(cljvec,3,numel(K2));
    
    Ri=rotations(:,:,k1);
    cli=Ri*cij; % R_{i}c_{ij}
    e=sum(cli.*cljvec);
    e=e(:);
    K2=K2(:);
    clerr(idx+1:idx+numel(nnzidx),:)=[repmat(k1,numel(nnzidx),1) K2(nnzidx) e(nnzidx)];
    idx=idx+numel(nnzidx);
end
clerr=clerr(1:idx,:);

% Make sure that clerr(:,3) is in [-1,1] (roundoff errors are allowed).
cosvals=clerr(:,3);

if norm(imag(cosvals))~=0
        idx=find(imag(cosvals)~=0);
        warning('Found %d cosines with imaginary components', numel(idx));
end
cosvals=real(cosvals);
    
coserr=cosvals(abs(cosvals)>1);
if ~isempty(coserr)
    if max(coserr)>1.0e-13
        idx=find(abs(cosvals)>1);
        warning('Found %d cosines outside [-1,1]', numel(idx));
    end;
end
cosvals=max(cosvals,-1);
cosvals=min(cosvals,1);
clerr(:,3)=real(cosvals);
clerr(:,3)=acosd(clerr(:,3)); % Convert to degrees


%% This is the old reference code before optimiation
% K=size(clmatrix,1);
% clerr=zeros(K*(K-1)/2,3); % Common lines embedding errors for each of the
% % K*(K-1)/2 common lines pairs. Each row of this array contains
% % (i,j,e), where i and j are the indices of the common lines, and e is
% % the dot product between R_{i}c_{ij} and R_{j}c_{ji}.
% idx=0;
% for k1=1:K-1
%     %fprintf('Processing k1=%d\n',k1);
%     for k2=k1+1:K
%         
%         if clmatrix(k1,k2)~=0 % Process only common lines that have not been
%             %rejected in previous iterations (if there were any such
%             %iterations), or during voting.
%             alphai=2*pi*(clmatrix(k1,k2)-1)/L;
%             alphaj=2*pi*(clmatrix(k2,k1)-1)/L;
%             cij=[cos(alphai) sin(alphai) 0].';
%             cji=[cos(alphaj) sin(alphaj) 0].';
%             cli=rotations(:,:,k1)*cij; % R_{i}c_{ij}
%             clj=rotations(:,:,k2)*cji; % R_{j}c_{ji}.
%             e=sum(cli.*clj);
%             idx=idx+1;
%             clerr(idx,:)=[k1 k2 e];
%         end
%     end
% end
% clerr=clerr(1:idx,:);
% clerr(:,3)=acosd(clerr(:,3)); % Convert to degrees
