function [idx1,idx2,Z3]=commonline_q(q,k1,k2,n_theta)
%
% Given a Kx4 array of quaternions, find the common lines between
% projections k1 and k2. Each projection has n_theta Fourier rays.
%
% idx1 is the index of the common line in projection k1. idx2 is the index
% of the common line in projection k2. Z3 is the direction vector of the
% common line.
%
% Yoel Shkolnisky, October 2008.

R1=q_to_rot(q(:,k1));
R2=q_to_rot(q(:,k2));

% Convert from Euler angles to a 3x3 rotation matrix.
R1=inv(R1);
R2=inv(R2);

% Rotated coordinate system is the columns of the rotation matrix.
% We use explicit multiplication to show which column corresponds to x,
% which to y, and which to z. See commonline_euler for the difference.
X1=R1*([1 0 0].');
Y1=R1*([0 1 0].');
Z1=R1*([0 0 1].');

X2=R2*([1 0 0].');
Y2=R2*([0 1 0].');
Z2=R2*([0 0 1].');

Z3=[Z1(2)*Z2(3)-Z1(3)*Z2(2);...
    Z1(3)*Z2(1)-Z1(1)*Z2(3);...
    Z1(1)*Z2(2)-Z1(2)*Z2(1)];

% Make sure the projections are not too close.
if norm(Z3)<1.0e-8
    warning('GCAR:normTooSmall','Images have same orientation');
end

Z3=Z3./norm(Z3);

% Compute coordinates of the common-line in each local coordinate system.
XY1=[X1 Y1];
XY2=[X2 Y2];
c1=(Z3.')*XY1;
c2=(Z3.')*XY2;

% Verify that the common-line is indeed common to both planes. The
% following warning should never happen! Just to make sure nothing went
% terribly wrong.  
ev1=XY1*c1(:)-Z3;
ev2=XY2*c2(:)-Z3;

if (norm(ev1)/norm(Z3)>1.0e-12) || (norm(ev2)/norm(Z3)>1.0e-12)
    warning('GCAR:largeErrors',...
        'Common line is not common. Error1 = %e, Error2 = %e',...
        norm(ev1)/norm(Z3),norm(ev2)/norm(Z3));
end

% Compute angle of the common line at each projection's coordinate system
theta1=atan2(c1(2),c1(1));
theta2=atan2(c2(2),c2(1));

PI=4*atan(1.0);
theta1=theta1+PI; % Shift from [-pi,pi] to [0,2*pi].
theta2=theta2+PI;

idx1=theta1/(2*PI)*n_theta;
idx2=theta2/(2*PI)*n_theta;
 
idx1=mod(round(idx1),n_theta);
idx2=mod(round(idx2),n_theta);