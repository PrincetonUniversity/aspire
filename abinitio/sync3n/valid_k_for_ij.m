
function phis = valid_k_for_ij ( clmatrix , i , j , L )
%VALID_k_FOR_ij For a given pair i,j, find all indeces k that are valid
% to be used for angle evaluation between i,j. Validity is determined
% according to singularity of the calculations based on common lines.
%
% Angle between i and j induced by each third projection k.
% 
% output:
% phis - n X 2 array of the cosine of the angle between i,j according to every k;
% and the index k the projection that creates that angle.
% 

% Configuration
TOL_idx = 1e-12;

% Initialization
N = size(clmatrix,1);
K = 1:N;
%K = K( K~=i & K~=j );
K = K( clmatrix(i,K)~=0 & clmatrix(j,K)~=0 ); % i,i & j,j assumed to be 0. other bad entries might be 0 as well.

% thetai = angle on Ci created by its intersection with Cj and Ck.
theta1 = (clmatrix(i,K)-clmatrix(i,j)) * 2*pi/L;
theta2 = (clmatrix(j,K)-clmatrix(j,i)) * 2*pi/L;
theta3 = (clmatrix(K,j)-clmatrix(K,i))' * 2*pi/L;

% cosines
c1 = cos(theta1);
c2 = cos(theta2);
c3 = cos(theta3);

% Each common-line corresponds to a point on the unit sphere. Denote the
% coordinates of these points by (Pix, Piy Piz), and put them in the matrix
%   M=[ P1x  P2x  P3x ; ...
%       P1y  P2y  P3y ; ...
%       P1z  P2z  P3z ].
%
% Then the matrix
%   C=[ 1 c1 c2 ;...
%       c1 1 c3 ;...
%       c2 c3 1],
% where c1,c2,c3 are given above, is given by C=M.'*M.
% For the points P1,P2, and P3 to form a triangle on the unit shpere, a
% necessary and sufficient condition is for C to be positive definite. This
% is equivalent to
%      1+2*c1*c2*c3-(c1^2+c2^2+c3^2)>0.
% However, this may result in a traingle that is too flat, that is, the
% angle between the projections is very close to zero. We therefore use the
% condition below
%       1+2*c1*c2*c3-(c1^2+c2^2+c3^2) > 1.0e-5
% This ensures that the smallest singular value (which is actually
% controlled by the determinant of C) is big enough, so the matrix is far
% from singular. This condition is equivalent to computing the singular
% values of C, followed by checking that the smallest one is big enough.
valid_triangles = find( 1+2*c1.*c2.*c3-(c1.^2+c2.^2+c3.^2) > 1e-5 );
theta1 = theta1(valid_triangles);
theta2 = theta2(valid_triangles);
c1 = c1(valid_triangles);
c2 = c2(valid_triangles);
c3 = c3(valid_triangles);

% Calculate Angles Cosines
cos_phi2 = (c3-c1.*c2)./(sqrt(1-c1.^2).*sqrt(1-c2.^2));

% Some synchronization must be applied when common line is
% out by 180 degrees.
% Here fix the angles between c_ij(c_ji) and c_ik(c_jk) to be smaller than pi/2,
% otherwise there will be an ambiguity between alpha and pi-alpha.
ind1 = ( (theta1 > pi+TOL_idx) | (theta1<-TOL_idx & theta1>-pi) );
ind2 = ( (theta2 > pi+TOL_idx) | (theta2<-TOL_idx & theta2>-pi) );
align180 = (ind1 & ~ind2 ) | (~ind1 & ind2 );
cos_phi2(align180) = -cos_phi2(align180);

% Fix numerical problems
bad_numeric = find(abs(cos_phi2) > 1);
if numel(find( (abs(cos_phi2(bad_numeric))-1) > 1.0e-12 )) > 0
    warning('GCAR:numericalProblem','cos_phi2>1 by too much... Setting to 1');
end
cos_phi2(bad_numeric) = sign(cos_phi2(bad_numeric));

% Return PHIS - cosines and good k indices
phis = [cos_phi2' K(valid_triangles)'];

end
