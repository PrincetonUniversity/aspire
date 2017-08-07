function [l_ij,l_ji]=commonline_R(Ri,Rj,L)
%
% Compute the common line induced by rotation matrixces Ri and Rj.
%
% Returns the indices of the common lines between images  (rotations) i and
% j in image i (l_ij) and in image j (l_ji), respectively.
% 
% NOTE that this function returns indices l_ij and l_ji that are
% zero-based, while commonlines matrices are 1 based (with 0 indicatiing
% that a common line has not been detected). To convert the returned
% indices l_ij and l_ji into 1-based, for example, to match the output of
% clmatrix_cheat, you must explicitly set
%   l_ij=l_ij+1; l_ji=l_ji+1; 
%
% See commonline_q for more information.
%
% Yoel Shkolnisky, July 2013.
%
% Revisions:
% Y.S. July 28, 2015    Revise header comment.

Ut=Rj*Ri.';

alphaij=atan2(Ut(3,1),-Ut(3,2));
alphaji=atan2(-Ut(1,3),Ut(2,3));

PI=4*atan(1.0);
alphaij=alphaij+PI; % Shift from [-pi,pi] to [0,2*pi].
alphaji=alphaji+PI;

l_ij=alphaij/(2*PI)*L;
l_ji=alphaji/(2*PI)*L;
 
l_ij=mod(round(l_ij),L);
l_ji=mod(round(l_ji),L);

%% Reference code (slower version of the above code):
% Ri=Ri.';
% Rj=Rj.';
% U=Ri.'*Rj;
% 
% alphaij=atan2(U(1,3),-U(2,3));
% alphaji=atan2(-U(3,1),U(3,2));
% 
% PI=4*atan(1.0);
% alphaij=alphaij+PI; % Shift from [-pi,pi] to [0,2*pi].
% alphaji=alphaji+PI;
% 
% l_ij=alphaij/(2*PI)*L;
% l_ji=alphaji/(2*PI)*L;
%  
% l_ij=mod(round(l_ij),L);
% l_ji=mod(round(l_ji),L);
