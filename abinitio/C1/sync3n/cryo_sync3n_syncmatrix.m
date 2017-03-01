function S = cryo_sync3n_syncmatrix(Rij)
%SYNC_MATRIX build 3NX3N synchronization matrix (without blocks weights).
% 
% Input:
% Rij 	For every pair i,j - the 3X3 matrix Ri'*Rj.
% 
% Output:
% S     3NX3N synchronization matrix (without blocks weights) where the
%       (i,j)'th 3x3 block is Ri'*Rj.
%
% This function was renamed from
%       sync_matrix ( N, Rij )
%
% Yoel Shkolnisky, March 2016.
% based on code by Ido Greenberg, 2015.

N=(1+sqrt(1+8*size(Rij,3)))/2; % Extract N from the number of relative rotaitons in Rij.
    % The number of relative rotation for a given N is (N choose 2) so,
    % solve for N from this number.
assert(N==round(N));  % Make sure we got an integer.

S = zeros(3*N);
I = eye(3);
idx = 0; % pair index
for i = 1:N  
    S( (3*i-2):(3*i) , (3*i-2):(3*i) ) = I; % Rii = I    
    for j = (i+1):N        
        idx = idx + 1;
        S( (3*i-2):(3*i) , (3*j-2):(3*j) ) = Rij(:,:,idx); % Rij
        S( (3*j-2):(3*j) , (3*i-2):(3*i) ) = Rij(:,:,idx)'; % Rji = Rij'        
    end
end
end
