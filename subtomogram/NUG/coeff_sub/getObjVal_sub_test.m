%{
INPUT:
    w: weights for numerical integration
    P: projections in Fourier space
    theta: angle discretization
    alpha: 1st Euler angle
    gamma: 3rd Euler angle
OUTPUT:
    objVal: value of the objective function at (alpha,gamma)X(i,j)

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function objVal = getObjVal_sub_test( w, P , alpha , beta, gamma , K_sub , ns, b_sub, n_theta, numCompThreads )

% [ref_clstack,~]=clmatrix_cheat_q(ref_q,n_theta);
% R_true = zeros(3,3,K_sub);
% for i = 1:K_sub
%     R_true(:,:,i) = q_to_rot(subq(:,i));
% end
% 
% Rij = R_true(:,:,2)*R_true(:,:,3)';
% [gij,bij,aij] = rot_to_EulerZYZ(Rij); % alpha and gamma are switched in Yutong's function
% k = 0;
% l = 0;
% bl = -pi/24*l;
% bk = -pi/24*k;
% s = sin([aij, bij, gij]);
% c = cos([aij, bij, gij]);
% sl = sin(bl);
% cl = cos(bl);
% sk = sin(bk);
% ck = cos(bk);
% 
% theta_ij = atan2((-c(3)*s(1)*sl - c(1)*c(2)*s(3)*sl + c(1)*s(2)*cl), ...
%     (-ck*(c(1)*c(3)*sl - c(2)*s(1)*s(3)*sl + s(1)*s(2)*cl) + sk*(s(2)*s(3)*sl+c(2)*cl)));
% theta_ji = atan2((-c(1)*s(3)*sk - c(2)*c(3)*s(1)*sk + c(3)*s(2)*ck), ...
%     (cl*(c(1)*c(3)*sk - c(2)*s(1)*s(3)*sk + s(2)*s(3)*ck) - sl*(s(1)*s(2)*sk+c(2)*ck))); 
% c_ij = round(mod(theta_ij,pi)/(2*pi/n_theta))+1;
% c_ji = round(mod(theta_ji,pi)/(2*pi/n_theta))+1;
% c_ij_true = round(mod(aij-pi/2,pi)/(2*pi/n_theta))+1;
% c_ji_true = round(mod(-gij-pi/2,pi)/(2*pi/n_theta))+1;
% 
% 
% 
% % theta_i = gamma - pi/2 and theta_j = -alpha - pi/2......THIS HAS BEEN VERIFIED
% [ ~ , idx_I ] = min( transpose( abs( repmat(theta,length(gamma),1) - mod( repmat(gamma,1,length(theta)) - pi/2 , 2*pi ) ) ) );
% [ ~ , idx_J ] = min( transpose( abs( repmat(theta,length(alpha),1) - mod( -repmat(alpha,1,length(theta)) - pi/2 , 2*pi ) ) ) );

% build objVal sub
count = 0;
objVal_subI = zeros( 1 , round(K_sub*(K_sub-1)/2) );
objVal_subJ = zeros( 1 , round(K_sub*(K_sub-1)/2) );
for j = 1:(K_sub-1)
for i = (j+1):K_sub
    count = count + 1;
    objVal_subI(count) = i;
    objVal_subJ(count) = j;
end
end
objVal_tiltI = repmat(0:ns-1,ns,1);
objVal_tiltJ = repmat((0:ns-1)',1,ns);
objVal_tiltI = objVal_tiltI(:);
objVal_tiltJ = objVal_tiltJ(:);

% get objVal entries corresponding to sub
alpha = alpha(:);
beta = beta(:);
gamma = gamma(:);
siz_so3 = length(alpha);
objVal = zeros( siz_so3 , length(objVal_subI) );
P = reshape(P,[],n_theta*ns,K_sub);

%parpool('local',numCompThreads);
for m = 1:siz_so3
    
    %display(m);
    [clmatrix_j, clmatrix_i] = getCL_sub_test( b_sub, [alpha(m), beta(m), gamma(m)], n_theta, ns ); %10/3
    
    clmatrix_i = clmatrix_i'; % so that it would be with objVal_tiltI;
    objVal_ij = abs( P(:, clmatrix_i(:)+objVal_tiltI*n_theta, objVal_subI) ...
        -  P(:, clmatrix_j(:)+objVal_tiltJ*n_theta, objVal_subJ));
    objVal_ij_sub = squeeze(sum( objVal_ij, 2)); % sum over tilts
    objVal_ij_rad = squeeze(sum( objVal_ij_sub.* repmat( w , 1 , size(objVal_ij,3) )  , 1 )); % sum over radius
    
    objVal( m , : ) = squeeze(objVal_ij_rad);
    
end
%delete(gcp);
% for i = 1:ns
%     tiltq(:,i+ns) = rot_to_q(q_to_rot(tiltq(:,i))*EulerZYZ_to_rot(alpha(m),beta(m),gamma(m))');
% end
% [ref_clstack,~]=clmatrix_cheat_q(tiltq,n_theta);

% objvalij = 0;
% for k = 1:ns
%     for l = 1:ns
%         cij = clmatrix_i(l,k);
%         cji = clmatrix_j(l,k);
%         objvalij = objvalij+sum(abs(P(:,cij+(k-1)*n_theta/2,2)-P(:,cji+(l-1)*n_theta/2,1)));
%     end
% end
% 
% [g12,b12,a12] = rot_to_EulerZYZ(Rij');
% mM = 1000;
% for m = 1:siz_so3;
%     dis = norm([alpha(m),beta(m),gamma(m)]-[a12,b12,g12]);
%     if dis < mM
%         mM = dis;
%         mI = m;
%     end
% end

end









