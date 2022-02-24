function [gR, scl_inds] = cryo_TO_group_elements(symmetry)

% Defines the symmetry group elements. 
%
% Note: there is a certain order of the symmetry group elements, which is 
% important for computing the self common lines set.
%
%   Input:
%       symmetry -      either 'T' or 'O'.
%
%   Output:
%       gR -            symmetry group elements: 3 x 3 x group order.    
%       n_gR -          group order, also the number of common lines.
%       n_scl_pairs -   number of self common lines.
%       
%   Written by Adi Shasha January 2021. 

% Set group elements for 'T':
if symmetry == 'T'
    n_gR = 12;                                  % The size of group T is 12.
    scl_inds = [2,4,6,8,10,11,12];
    
    gR = zeros(3,3,n_gR);                       % T symmetry group.
    gR(:,:,01)  = eye(3);
    gR(:,:,02) = [ 0  0  1; 1  0  0; 0  1  0];      % axis: [ 1, 1, 1] angle: 120
    gR(:,:,03) = [ 0  1  0; 0  0  1; 1  0  0];      % axis: [ 1, 1, 1] angle: 240
    gR(:,:,04) = [ 0  0 -1; 1  0  0; 0 -1  0];      % axis: [-1,-1, 1] angle: 120
    gR(:,:,05) = [ 0  1  0; 0  0 -1;-1  0  0];      % axis: [-1,-1, 1] angle: 240
    gR(:,:,06) = [ 0  0 -1;-1  0  0; 0  1  0];      % axis: [ 1,-1,-1] angle: 120
    gR(:,:,07) = [ 0 -1  0; 0  0  1;-1  0  0];      % axis: [ 1,-1,-1] angle: 240
    gR(:,:,08) = [ 0  0  1;-1  0  0; 0 -1  0];      % axis: [-1, 1,-1] angle: 120
    gR(:,:,09) = [ 0 -1  0; 0  0 -1; 1  0  0];      % axis: [-1, 1,-1] angle: 240
    gR(:,:,10) = [ 1  0  0; 0 -1  0; 0  0 -1];      % axis: [ 1, 0, 0] angle: 180
    gR(:,:,11) = [-1  0  0; 0  1  0; 0  0 -1];      % axis: [ 0, 1, 0] angle: 180
    gR(:,:,12) = [-1  0  0; 0 -1  0; 0  0  1];      % axis: [ 0, 0, 1] angle: 180

% Set group elements for 'O':
elseif symmetry == 'O'
    n_gR  = 24;                                 % The size of group O is 24.
    scl_inds = [2,4,6,8,10,11,12,14,16,17,18,19,20,21,22,23,24];
    
    gR = zeros(3,3,n_gR);                       % O symmetry group.
    gR(:,:,01) = [ 0 -1  0; 1  0  0; 0  0  1];   
    gR(:,:,02) = [ 0  1  0;-1  0  0; 0  0  1];      % gR_1^T = gR_2
    gR(:,:,03) = [ 1  0  0; 0  0 -1; 0  1  0];   
    gR(:,:,04) = [ 1  0  0; 0  0  1; 0 -1  0];      % gR_3^T = gR_4
    gR(:,:,05) = [ 0 -1  0; 0  0 -1; 1  0  0];   
    gR(:,:,06) = [ 0  0  1;-1  0  0; 0 -1  0];      % gR_5^T = gR_6
    gR(:,:,07) = [ 0 -1  0; 0  0  1;-1  0  0];   
    gR(:,:,08) = [ 0  0 -1;-1  0  0; 0  1  0];      % gR_7^T = gR_8
    gR(:,:,09) = [ 0  1  0; 0  0  1; 1  0  0];   
    gR(:,:,10) = [ 0  0  1; 1  0  0; 0  1  0];      % gR_9^T = gR_10
    gR(:,:,11) = [ 0  0  1; 0  1  0;-1  0  0];   
    gR(:,:,12) = [ 0  0 -1; 0  1  0; 1  0  0];      % gR_11^T = gR_12
    gR(:,:,13) = [ 0  1  0; 0  0 -1;-1  0  0];   
    gR(:,:,14) = [ 0  0 -1; 1  0  0; 0 -1  0];      % gR_13^T = gR_14
    gR(:,:,15) = [ 1  0  0; 0  1  0; 0  0  1];   
    gR(:,:,16) = [-1  0  0; 0 -1  0; 0  0  1]; 
    gR(:,:,17) = [-1  0  0; 0  0 -1; 0 -1  0];   
    gR(:,:,18) = [ 1  0  0; 0 -1  0; 0  0 -1];
    gR(:,:,19) = [ 0 -1  0;-1  0  0; 0  0 -1];   
    gR(:,:,20) = [-1  0  0; 0  1  0; 0  0 -1];
    gR(:,:,21) = [ 0  1  0; 1  0  0; 0  0 -1];   
    gR(:,:,22) = [-1  0  0; 0  0  1; 0  1  0];
    gR(:,:,23) = [ 0  0  1; 0 -1  0; 1  0  0];   
    gR(:,:,24) = [ 0  0 -1; 0 -1  0;-1  0  0];

else
   error('Wrong symmetry type: should be T or O.'); 
end