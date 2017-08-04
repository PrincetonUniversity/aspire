function test_est_orientations(K,L,p)

% K is the number of projections
% L is the number of radial lines within a projection
% p is the probability that correlating two projection images yields the
%    correct common line



%%%%%% Generate K random rotations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                     %%
% The 3-sphere S^3 in R^4 is a double cover of the rotation group SO(3)
% SO(3) = RP^3
% We identify unit norm quaternions a^2+b^2+c^2+d^2=1 with group elements
% The antipodal points (-a,-b,-c,-d) and (a,b,c,d) are identified as the
% same group elements

rot_matrices = rand_rots(K);

% calculate inverse rotation matrices (just transpose)
inv_rot_matrices = zeros(3,3,K);

inv_rot_matrix = zeros(3);
for k=1:K;
    rot_matrix = rot_matrices(:,:,k);
    inv_rot_matrix = rot_matrix'; % inv(R)=R^T
    inv_rot_matrices(:,:,k) = inv_rot_matrix;
    inv_rot_matrix(:,:,k) = inv(rot_matrices(:,:,k));
end;

n_x(:,:) = inv_rot_matrices(:,1,:);
n_y(:,:) = inv_rot_matrices(:,2,:);
n_z(:,:) = inv_rot_matrices(:,3,:);

fprintf('Computing common lines... ');
% Find common lines -- the intersection points of the circles
common_lines_matrix = zeros(K);

eqs_matrix = zeros(3,4);

for k1=1:K;
    for k2=(k1+1):K;
        
        %%%%%%%%%%%%% find common lines using linear algebra %%%%%%%%%%%%%
        %
        % find the common line of two circles
        % plane(k1) = a1 * plane1(k1) + b1 * plane2(k1)
        % plane(k2) = a2 * plane1(k2) + b2 * plane2(k2)
        % common line: plane(k1) = plane(k2) ==>
        % a1 * plane1(k1) + b1 * plane2(k1) = a2 * plane1(k2) + b2 * plane2(k2)
        % or equivaletly
        % plane1_x(k1)*a1+plane2_x(k1)*b1-plane1_x(k2)*a2-plane2_x(k2)*b2 = 0
        % plane1_y(k1)*a1+plane2_y(k1)*b1-plane1_y(k2)*a2-plane2_y(k2)*b2 = 0
        % plane1_z(k1)*a1+plane2_z(k1)*b1-plane1_z(k2)*a2-plane2_z(k2)*b2 = 0
        % 3 equations in 4 unknowns a1,b1,a2,b2.
        eqs_matrix = [n_x(1,k1), n_y(1,k1), -n_x(1,k2), -n_y(1,k2) ; ...
            n_x(2,k1), n_y(2,k1), -n_x(2,k2), -n_y(2,k2) ; ...
            n_x(3,k1), n_y(3,k1), -n_x(3,k2), -n_y(3,k2)];
        Z = 2*null(eqs_matrix); % the null space is the special solution === direction of common line
        % the factors 2 is because matlab returns Z normalized, that is, |Z|^2=1
        % and we want a1^2+b1^2 = a2^2 + b2^2 = 1
        
        % now, we know a1=cos(phi(k1))and b1=sin(phi(k1))
        % so, phi(k1) = atan(b1/a1)
        phi_1 = atan2(Z(2),Z(1));
        if (phi_1 < 0)
            phi_1 = phi_1 + 2*pi;
        end;
        % and phi(k2) = atan(b2/a2)
        phi_2 = atan2(Z(4),Z(3));
        if (phi_2 < 0)
            phi_2 = phi_2 + 2*pi;
        end;
        l1 = mod(round(phi_1*L/(2*pi))+L-1,L)+1;
        l2 = mod(round(phi_2*L/(2*pi))+L-1,L)+1;
        
        common_lines_matrix(k1,k2) = l1;
        common_lines_matrix(k2,k1) = l2;
    end;
end;

fprintf('Finished!\n');

%%%%% Perturb the common lines matrix
%%
is_perturbed = zeros(K); % to verify the success of the algorithm we store which common lines were perturbed
for k1=1:K;
    for k2=(k1+1):K;
        r = rand;
        % with probability 1-p the common line needs to be perturbed
        if (r > p)
            is_perturbed(k1,k2) = 1;
            common_lines_matrix(k1,k2) = floor(rand*L)+1;
            common_lines_matrix(k2,k1) = floor(rand*L)+1;
        end;
    end;
end;

%%%%% Test different orientation determination algorithms
MSEs = zeros(6,1);
Time = zeros(6,1);
%% Test est_orientations_LS
tic;
est_inv_rots = est_orientations_LS(common_lines_matrix, L);
Time(1) = toc;
MSEs(1) = check_MSE(est_inv_rots, rot_matrices);

% With spectral norm constraint
pars.alpha = 2/3;
tic;
est_inv_rots = est_orientations_LS(common_lines_matrix, L, pars);
Time(2) = toc;
MSEs(2) = check_MSE(est_inv_rots, rot_matrices);

%% Test est_orientations_LUD
% ADMM
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,L);
Time(3) = toc;
MSEs(3) = check_MSE(est_inv_rots, rot_matrices);

% ADMM with spectral norm constraint
pars.alpha = 2/3;
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,L, pars);
Time(4) = toc;
MSEs(4) = check_MSE(est_inv_rots, rot_matrices);

% IRLS
pars.solver = 'IRLS';
pars.alpha = 0;
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,L, pars);
Time(5) = toc;
MSEs(5) = check_MSE(est_inv_rots, rot_matrices);

% IRLS with spectral norm constraint
pars.alpha = 2/3;
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,L, pars);
Time(6) = toc;
MSEs(6) = check_MSE(est_inv_rots, rot_matrices);

%% Print the MSEs and cost time of the results
fprintf('K = %d, L = %d, p = %f\n', K, L, p);
fprintf('=======================================================================================\n')
fprintf(' Exp    LS     LUD     alpha     SDPLR     ADMM      IRLS     MSE             Time\n')
fprintf('  1     Y                          Y                          %1.5f        %6.2f\n', MSEs(1), Time(1));
fprintf('  2     Y               2/3                 Y                 %1.5f        %6.2f\n', MSEs(2), Time(2));
fprintf('  3             Y                           Y                 %1.5f        %6.2f\n', MSEs(3), Time(3));
fprintf('  4             Y       2/3                 Y                 %1.5f        %6.2f\n', MSEs(4), Time(4));
fprintf('  5             Y                                     Y       %1.5f        %6.2f\n', MSEs(5), Time(5));
fprintf('  6             Y       2/3                           Y       %1.5f        %6.2f\n', MSEs(6), Time(6));
fprintf('=======================================================================================\n')

% For example:
% K = 500, L = 360, p = 0.500000
% =======================================================================================
%  Exp    LS     LUD     alpha     SDPLR     ADMM      IRLS     MSE             Time
%   1     Y                          Y                          0.01430          2.26
%   2     Y               2/3                 Y                 0.01584         24.07
%   3             Y                           Y                 0.00150        109.29
%   4             Y       2/3                 Y                 0.00395        138.04
%   5             Y                                     Y       0.00061         32.16
%   6             Y       2/3                           Y       0.00274        1075.34
% =======================================================================================
