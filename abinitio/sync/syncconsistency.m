function clerr=syncconsistency(rotations,clmatrix,L)
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


K=size(clmatrix,1);
clerr=zeros(K*(K-1)/2,3); % Common lines embedding errors for each of the
% K*(K-1)/2 common lines pairs. Each row of this array contains
% (i,j,e), where i and j are the indices of the common lines, and e is
% the dot product between R_{i}c_{ij} and R_{j}c_{ji}.
idx=0;
for k1=1:K-1
    %fprintf('Processing k1=%d\n',k1);
    for k2=k1+1:K
        
        if clmatrix(k1,k2)~=0 % Process only common lines that have not been
            %rejected in previous iterations (if there were any such
            %iterations), or during voting.
            alphai=2*pi*(clmatrix(k1,k2)-1)/L;
            alphaj=2*pi*(clmatrix(k2,k1)-1)/L;
            cij=[cos(alphai) sin(alphai) 0].';
            cji=[cos(alphaj) sin(alphaj) 0].';
            cli=rotations(:,:,k1)*cij; % R_{i}c_{ij}
            clj=rotations(:,:,k2)*cji; % R_{j}c_{ji}.
            e=sum(cli.*clj);
            idx=idx+1;
            clerr(idx,:)=[k1 k2 e];
        end
    end
end
clerr=clerr(1:idx,:);
clerr(:,3)=acosd(clerr(:,3)); % Convert to degrees
