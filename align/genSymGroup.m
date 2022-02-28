function G = genSymGroup(sym)
% This function generates the symmetry group of the given symmetry type. 

% The symmetry elements are being generated according to the most common 
% coordinate systems of molecules from the EMDB. Note that it is necessary  
% to check that the generated symmetry group is indeed the appropriate one.  
%% Input:
% sym- the symmetry type- 'Cn'\'Dn'\'T'\'O'\'I', where n is the the symmetry
%      order.  
%% Output:
% G- Array of matrices of size 3x3xn containing the symmetry group elemnts. 

%%

if ~exist('sym','var') || isempty(sym)
    error('sym must be an input')
end

symgroup = sym(1);
if numel(sym) > 1
        n_s = str2double(sym(2:end));
else
    n_s=0;
end

if symgroup == 'C' 
    if n_s == 0
        error('The symmetry order must be an input in case of cyclic symmetry')
    end
    G = zeros(3,3,n_s);                      % C symmetry group.
    theta = 2*180/n_s;
    for i = 1:n_s
        G(:,:,i) = rotz((i-1)*theta);
    end
    
elseif symgroup == 'D'
    if n_s == 0
        error('The symmetry order must be an input in case of dihedral symmetry')
    end
    G = zeros(3,3,2*n_s);                    % D symmetry group.
    theta = 2*180/n_s;
    for i = 1:n_s
        G(:,:,i) = rotz((i-1)*theta);
    end
    G(:,:,n_s+1) = roty(180);
    for i = 2:n_s
        G(:,:,n_s+i) = G(:,:,n_s+1)*G(:,:,i);
    end
    
elseif symgroup == 'T'
    n = 12;                                  % The size of group T is 12.
    G = zeros(3,3,n);                        % T symmetry group.
    G(:,:,1)  = eye(3);
    G(:,:,02) = [ 0  0  1; 1  0  0; 0  1  0];      % axis: [ 1, 1, 1] angle: 120
    G(:,:,03) = [ 0  1  0; 0  0  1; 1  0  0];      % axis: [ 1, 1, 1] angle: 240    G(03) = G(02)^T
    G(:,:,04) = [ 0  0 -1; 1  0  0; 0 -1  0];      % axis: [-1,-1, 1] angle: 120
    G(:,:,05) = [ 0  1  0; 0  0 -1;-1  0  0];      % axis: [-1,-1, 1] angle: 240    G(05) = G(04)^T
    G(:,:,06) = [ 0  0 -1;-1  0  0; 0  1  0];      % axis: [ 1,-1,-1] angle: 120
    G(:,:,07) = [ 0 -1  0; 0  0  1;-1  0  0];      % axis: [ 1,-1,-1] angle: 240    G(07) = G(06)^T
    G(:,:,08) = [ 0  0  1;-1  0  0; 0 -1  0];      % axis: [-1, 1,-1] angle: 120
    G(:,:,09) = [ 0 -1  0; 0  0 -1; 1  0  0];      % axis: [-1, 1,-1] angle: 240    G(09) = G(08)^T

    G(:,:,10) = [ 1  0  0; 0 -1  0; 0  0 -1];      % axis: [ 1, 0, 0] angle: 180    G(10) = G(10)^T
    G(:,:,11) = [-1  0  0; 0  1  0; 0  0 -1];      % axis: [ 0, 1, 0] angle: 180    G(11) = G(11)^T
    G(:,:,12) = [-1  0  0; 0 -1  0; 0  0  1];      % axis: [ 0, 0, 1] angle: 180    G(12) = G(12)^T
    
elseif symgroup == 'O'
    n  = 24;                                  % The size of group O is 24.
    G = zeros(3,3,n);                         % O symmetry group.
    G(:,:,01) = eye(3); 
    G(:,:,02) = [ 0 -1  0; 1  0  0; 0  0  1];   
    G(:,:,03) = [ 0  1  0;-1  0  0; 0  0  1];      % G_3 = G_2.'
    G(:,:,04) = [ 1  0  0; 0  0 -1; 0  1  0];   
    G(:,:,05) = [ 1  0  0; 0  0  1; 0 -1  0];      % G_5 = G_4.'
    G(:,:,06) = [ 0 -1  0; 0  0 -1; 1  0  0];   
    G(:,:,07) = [ 0  0  1;-1  0  0; 0 -1  0];      % G_7 = G_6.'
    G(:,:,08) = [ 0 -1  0; 0  0  1;-1  0  0];   
    G(:,:,09) = [ 0  0 -1;-1  0  0; 0  1  0];      % G_9 = G_8.'
    G(:,:,10) = [ 0  1  0; 0  0  1; 1  0  0];   
    G(:,:,11) = [ 0  0  1; 1  0  0; 0  1  0];      % G_11 = G_10.'
    G(:,:,12) = [ 0  0  1; 0  1  0;-1  0  0];   
    G(:,:,13) = [ 0  0 -1; 0  1  0; 1  0  0];      % G_13 = G_12.'
    G(:,:,14) = [ 0  1  0; 0  0 -1;-1  0  0];   
    G(:,:,15) = [ 0  0 -1; 1  0  0; 0 -1  0];      % G_15 = G_14.'

    G(:,:,16) = [-1  0  0; 0 -1  0; 0  0  1];      % 2-fold
    G(:,:,17) = [-1  0  0; 0  0 -1; 0 -1  0];   
    G(:,:,18) = [ 1  0  0; 0 -1  0; 0  0 -1];
    G(:,:,19) = [ 0 -1  0;-1  0  0; 0  0 -1];   
    G(:,:,20) = [-1  0  0; 0  1  0; 0  0 -1];
    G(:,:,21) = [ 0  1  0; 1  0  0; 0  0 -1];   
    G(:,:,22) = [-1  0  0; 0  0  1; 0  1  0];
    G(:,:,23) = [ 0  0  1; 0 -1  0; 1  0  0];   
    G(:,:,24) = [ 0  0 -1; 0 -1  0;-1  0  0];
    
elseif symgroup == 'I'    
    n  = 60;                                 % The size of group I is 60.
    G = zeros(3,3,n);                        % I symmetry group.
    phi = (1 + sqrt(5))/2;
    G(:,:,01) = [ 1, 0, 0; 0, 1, 0; 0, 0, 1];               % Identity (D2)
    % 6 rotation axis joining the extreme opposite vertices, by angles 
    % 2pi/5, 4pi/5, 6pi/5 and 8pi/5
    G(:,:,02) = axang2rotm([0 phi 1 (2*pi)/5]);
    G(:,:,03) = axang2rotm([0 phi 1 (4*pi)/5]);
    G(:,:,04) = axang2rotm([0 phi 1 (6*pi)/5]);
    G(:,:,05) = axang2rotm([0 phi 1 (8*pi)/5]);

    G(:,:,06) = axang2rotm([0 phi -1 (2*pi)/5]);
    G(:,:,07) = axang2rotm([0 phi -1 (4*pi)/5]);
    G(:,:,08) = axang2rotm([0 phi -1 (6*pi)/5]);
    G(:,:,09) = axang2rotm([0 phi -1 (8*pi)/5]);

    G(:,:,10) = axang2rotm([1 0 phi (2*pi)/5]);
    G(:,:,11) = axang2rotm([1 0 phi (4*pi)/5]);
    G(:,:,12) = axang2rotm([1 0 phi (6*pi)/5]);
    G(:,:,13) = axang2rotm([1 0 phi (8*pi)/5]);

    G(:,:,14) = axang2rotm([-1 0 phi (2*pi)/5]);
    G(:,:,15) = axang2rotm([-1 0 phi (4*pi)/5]);
    G(:,:,16) = axang2rotm([-1 0 phi (6*pi)/5]);
    G(:,:,17) = axang2rotm([-1 0 phi (8*pi)/5]);

    G(:,:,18) = axang2rotm([phi -1 0 (2*pi)/5]);
    G(:,:,19) = axang2rotm([phi -1 0 (4*pi)/5]);
    G(:,:,20) = axang2rotm([phi -1 0 (6*pi)/5]);
    G(:,:,21) = axang2rotm([phi -1 0 (8*pi)/5]);

    G(:,:,22) = axang2rotm([phi 1 0 (2*pi)/5]);
    G(:,:,23) = axang2rotm([phi 1 0 (4*pi)/5]);
    G(:,:,24) = axang2rotm([phi 1 0 (6*pi)/5]);
    G(:,:,25) = axang2rotm([phi 1 0 (8*pi)/5]);

    % 10 rotation axis joining the centers of opposite faces, by angles 
    % 2pi/3 and 4pi/3.  
    G(:,:,26) = axang2rotm([0 -1 phi^2 (2*pi)/3]);
    G(:,:,27) = axang2rotm([0 -1 phi^2 (4*pi)/3]);

    G(:,:,28) = axang2rotm([0 1 phi^2 (2*pi)/3]);
    G(:,:,29) = axang2rotm([0 1 phi^2 (4*pi)/3]);

    G(:,:,30) = axang2rotm([-1 phi^2 0 (2*pi)/3]);
    G(:,:,31) = axang2rotm([-1 phi^2 0 (4*pi)/3]);

    G(:,:,32) = axang2rotm([1 phi^2 0 (2*pi)/3]);
    G(:,:,33) = axang2rotm([1 phi^2 0 (4*pi)/3]);

    G(:,:,34) = axang2rotm([phi^2 0 -1 (2*pi)/3]);
    G(:,:,35) = axang2rotm([phi^2 0 -1 (4*pi)/3]);

    G(:,:,36) = axang2rotm([phi^2 0 1 (2*pi)/3]);
    G(:,:,37) = axang2rotm([phi^2 0 1 (4*pi)/3]);

    G(:,:,38) = axang2rotm([1 -1 1 (2*pi)/3]);
    G(:,:,39) = axang2rotm([1 -1 1 (4*pi)/3]);

    G(:,:,40) = axang2rotm([-1 1 1 (2*pi)/3]);
    G(:,:,41) = axang2rotm([-1 1 1 (4*pi)/3]);

    G(:,:,42) = axang2rotm([1 1 -1 (2*pi)/3]);
    G(:,:,43) = axang2rotm([1 1 -1 (4*pi)/3]);

    G(:,:,44) = axang2rotm([1 1 1 (2*pi)/3]);
    G(:,:,45) = axang2rotm([1 1 1 (4*pi)/3]);

    % 15 rotation axis joining the midpoints of opposite edges, by angle pi.
    G(:,:,46) = axang2rotm([0 1 0 pi]);
    G(:,:,47) = axang2rotm([0 0 1 pi]);
    G(:,:,48) = axang2rotm([1 0 0 pi]);
    G(:,:,49) = axang2rotm([-1/phi 1 phi pi]);
    G(:,:,50) = axang2rotm([1/phi -1 phi pi]);
    G(:,:,51) = axang2rotm([1/phi 1 -phi pi]);
    G(:,:,52) = axang2rotm([1/phi 1 phi pi]);
    G(:,:,53) = axang2rotm([-1 phi 1/phi pi]);
    G(:,:,54) = axang2rotm([1 -phi 1/phi pi]);
    G(:,:,55) = axang2rotm([1 phi -1/phi pi]);
    G(:,:,56) = axang2rotm([1 phi 1/phi pi]);
    G(:,:,57) = axang2rotm([-phi 1/phi 1 pi]);
    G(:,:,58) = axang2rotm([phi -1/phi 1 pi]);
    G(:,:,59) = axang2rotm([phi 1/phi -1 pi]);
    G(:,:,60) = axang2rotm([phi 1/phi 1 pi]);
    
    O_g = [0 1 0; 1 0 0; 0 0 1]; %XXXXX
    for i = 1:60 %XXXXX
        G(:,:,i) = O_g*G(:,:,i)*O_g.'; %XXXXX
    end %XXXXX
       
else 
    error('sym was not entered properly') 
end

    
              
  