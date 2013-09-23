function R=rotratio(cl,L)
% R=rotratio(cl,L)
%
% Given a 3x3 common lines matrix, where the index of each common line is
% between 1 and L, compute the rotation that takes image 1 to image 2.
%
% Error codes:
%    -101 Triangle too small
%    -102 Condition number too large (larger than MAX_COND)
%
% Yoel Shkolnisky, August 2010

EPS=1.0e-13; % Should be 1.0e-13 after fixing XXX below.
MAX_COND=1000; % Largest allowed condition number for the system of equtions

% Find Q12, Q13, Q23
a=cos(2*pi*(cl(3,2)-cl(3,1))/L);
b=cos(2*pi*(cl(2,3)-cl(2,1))/L);
c=cos(2*pi*(cl(1,3)-cl(1,2))/L);

% Make sure that the triangle is not too small. This will happen if the
% common line between (say) cl(1,2) is close to cl(1,3). 
% To eliminate that, we require that det(G)=1+2abc-(a^2+b^2+c^2) is large
% enough.

if 1+2*a*b*c-(a^2+b^2+c^2)<1.0e-5
    R=-101;
    return;
end

G = [1 a b;...
     a 1 c;...
     b c 1];
 
if det(G)<0
    aaa=1;
end

Q=chol(G);


assert(norm(sum(Q.^2)-1)<EPS);
assert(norm(Q.'*Q-G)<EPS);

Q23=Q(:,1);
Q13=Q(:,2);
Q12=Q(:,3);

%
% Find the rotation matrix that takes:
%   1) the common line cl(1,2) to Q12;
%   2) the common line cl(1,3) to Q13;
%   3) their normal (0,0,1) to the normal of Q12 and Q13.
%
alpha12=2*pi*(cl(1,2)-1)/L;
alpha13=2*pi*(cl(1,3)-1)/L;
n1=mycross(Q12,Q13); n1=n1./norm(n1);
U1(:,1)=[cos(alpha12) ; sin(alpha12) ; 0];
U1(:,2)=[cos(alpha13) ; sin(alpha13) ; 0];
U1(:,3)=[0 ; 0; 1];

if det(U1)<0
    U1(:,3)=-U1(:,3); % Make sure we have a right-handed system)
end

if (cond(U1)>MAX_COND) %Find an explicit expression for the condition number
    fprintf('cond(U1)=%e\n',cond(U1));
    R=-102;
    return;
end

Q1=[Q12 Q13 n1];
T1=Q1*inv(U1);

if norm(T1.'*T1)-1>MAX_COND*EPS
    aaa=1;
end

assert(norm(T1.'*T1)-1<MAX_COND*EPS);
assert(abs(det(T1)-1)<MAX_COND*EPS);

[U,S,V]=svd(T1);
R1=U*V.';

% Verify that we got a rotation matrix
%assert(norm(R1.'*R1)-1<EPS);
if norm(R1.'*R1)-1>EPS
    aaa=1;
end

if  abs(det(R1)-1)>EPS
    aaa=1;
end
%assert(abs(det(R1)-1)<EPS || abs(det(R1)+1)<EPS);
assert(abs(det(R1)-1)<EPS);

%
% Find the rotation matrix that takes:
%   1) the common line cl(2,1) to Q12;
%   2) the common line cl(2,3) to Q23;
%   3) their normal (0,0,1) to the normal of Q12 and Q23.
%
alpha21=2*pi*(cl(2,1)-1)/L;
alpha23=2*pi*(cl(2,3)-1)/L;
n2=mycross(Q12,Q23); n2=n2./norm(n2);

U2(:,1)=[cos(alpha21) ; sin(alpha21) ; 0];
U2(:,2)=[cos(alpha23) ; sin(alpha23) ; 0];
U2(:,3)=[0 ; 0; 1];

if det(U2)<0
    U2(:,3)=-U2(:,3); % Make sure we have a right-handed system)
end

if (cond(U2)>MAX_COND) %Find an explicit expression for the condition number
    fprintf('cond(U1)=%e\n',cond(U2));
    R=-102;
    return;
end

Q2=[Q12 Q23 n2];
T2=Q2*inv(U2);

assert(norm(T2.'*T2)-1<MAX_COND*EPS);
assert(abs(det(T2)-1)<MAX_COND*EPS);

[U,S,V]=svd(T2);
R2=U*V.';

% Verify that we got a rotation matrix.
%assert(norm(R2.'*R2)-1<EPS);
if (norm(R2.'*R2)-1>EPS)
    aaa=1;
end
    
%assert(abs(det(R2)-1)<EPS || abs(det(R2)+1)<EPS);
assert(abs(det(R2)-1)<EPS);

% Return the ratio of the two rotation matrices.
R=R1.'*R2;