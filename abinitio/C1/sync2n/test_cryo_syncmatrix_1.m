% Test the function cryo_syncmatrix.
% 
% All messages should print OK, with errors of the order of machine
% precision.
%
% Yoel Shkolnisky, October 2010

function test_cryo_syncmatrix_1

initstate;
TOL=1.0e-14;
K=50;
rots_ref = rand_rots(K);
L=1E15;     % Use a large number of lines per image, so we don't have discretization errors.
cl=clmatrix_cheat(rots_ref,L);
S=cryo_syncmatrix(cl,L,rots_ref);
s=eig(S);
s=sort(s);

% There should be only three large eigenvalue, and all the others very
% close to zero.
e1=norm(s(1:end-3)/norm(S)); 
if e1>TOL
    fprintf('**** Rank of S is larger than three e1=%e\n',e1);
else
    fprintf('Rank test OK. e1=%e\n',e1);
end

% Check that the MSE of the recovered orientations is very small
[~,diff,mse]=cryo_syncrotations(S,rots_ref);
e2=max(diff(:));

if e2>TOL
    fprintf('**** Maximal rotation is error is large e2=%e\n',e2);
else
    fprintf('Maximal rotation is error is OK e2=%e\n',e2);
end

if mse>TOL
    fprintf('**** mse is large mse=%e\n',mse);
else
    fprintf('mse OK mse=%e\n',mse);
end

