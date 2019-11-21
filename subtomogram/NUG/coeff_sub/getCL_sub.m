function [ clmatrix_i, clmatrix_j ] = getCL_sub( b_sub, eul_ij, n_theta, ns )
% Compute the common line indices for a fixed set of orientation between
% subtomograms
% 
% Input
% b_sub: the increment in angle with tilt series discretization
% sub_q: orientations of subtomograms
% n_theta: angle discretization 
% K_sub: number of subtomograms
% ns: number of tilts in each subtomogram
%
% Output
% clmatrix: matrix with (i,j) represent c_ij, common line on the i-th
% projection from the j-th image.

%b_sub = pi/24;

% R_true = zeros(3,3,K_sub);
% for i = 1:K_sub
%     R_true(:,:,i) = q_to_rot(subq(:,i));
% end

clmatrix_i = zeros(ns);
clmatrix_j = zeros(ns);
mid_s = (ns+1)/2;

% i = 1;
% j = 2;
% Rij = q_to_rot(subq(:,i))*q_to_rot(subq(:,j))';
% [gij,bij,aij] = rot_to_EulerZYZ(Rij); % alpha and gamma are switched in Yutong's function
% eul_ij = [aij,bij,gij];

s = sin(eul_ij);
c = cos(eul_ij);
for k = 1:ns
    for l = 1:ns
        bl = -b_sub*(l-mid_s);
        bk = -b_sub*(k-mid_s);
        sl = sin(bl);
        cl = cos(bl);
        sk = sin(bk);
        ck = cos(bk);
        theta_ij = atan2((-c(3)*s(1)*sl - c(1)*c(2)*s(3)*sl + c(1)*s(2)*cl), ...
            (-ck*(c(1)*c(3)*sl - c(2)*s(1)*s(3)*sl + s(1)*s(2)*cl) + sk*(s(2)*s(3)*sl+c(2)*cl)));
        theta_ji = atan2((-c(1)*s(3)*sk - c(2)*c(3)*s(1)*sk + c(3)*s(2)*ck), ...
            (cl*(c(1)*c(3)*sk - c(2)*s(1)*s(3)*sk + s(2)*s(3)*ck) - sl*(s(1)*s(2)*sk+c(2)*ck))); 
        c_ij = round(mod(theta_ij,pi)/(2*pi/n_theta))+1;
        c_ji = round(mod(theta_ji,pi)/(2*pi/n_theta))+1;
        if c_ij > n_theta/2
            c_ij = mod(c_ij, n_theta/2);
        end
        if c_ji > n_theta/2
            c_ji = mod(c_ji, n_theta/2);
        end
        clmatrix_i(k,l) = n_theta/2-c_ij+2;
        clmatrix_j(l,k) = c_ji;
    end
end
%     end
% end
%clmatrix_j = clmatrix_j';

% [ref_clstack,~]=clmatrix_cheat_q(ref_q,n_theta);
% diff = mod(clmatrix-ref_clstack,n_theta/2); % checked with reference!

end

