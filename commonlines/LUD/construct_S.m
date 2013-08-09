function [S] = construct_S(common_lines_matrix, L)
% Compute the 2K x 2K commonline matrix S
%
% Input:
%   common_lines_matrix: a K x K matrix. (k1,k2) and (k2,k1) contain the index
%       of the common line of projections k1 and k2. (k1,k2)  contains the
%       index of the common line in the projection k1. (k2,k1) contains the
%       index of the common line in k2. 
%
%   L:  the number of radial lines within the Fourier transform of a
%       projection.
%
% Output: 
%   S:  2K x 2K commonline matrix.
%
% Lanhui Wang, Aug 8, 2013

K = size(common_lines_matrix,1); % K is the number of projections

S11 = zeros(K);
S12 = zeros(K);
S21 = zeros(K);
S22 = zeros(K);

for k1=1:K;
    k2=(k1+1):K;
    l1 = common_lines_matrix(k1,k2)-1;
    l2 = common_lines_matrix(k2,k1)-1;
    l1=l1(:);
    l2=l2(:);
    
    x12 = cos(2*pi*l1/L);
    y12 = sin(2*pi*l1/L);
    x21 = cos(2*pi*l2/L);
    y21 = sin(2*pi*l2/L);
    
    S11(k1,k2) = x12.*x21;
    S11(k2,k1) = x21.*x12;
    
    S12(k1,k2) = x12.*y21;
    S12(k2,k1) = x21.*y12;
    
    S21(k1,k2) = y12.*x21;
    S21(k2,k1) = y21.*x12;
    
    S22(k1,k2) = y12.*y21;
    S22(k2,k1) = y21.*y12;
    
    
end;

S = [S11 S12 ; S21 S22]; 