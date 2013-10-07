function ms=CircSymmetrize(m)
% function ms=CircSymmetrize(m)
% Given an n x n image m, compute its rotational average about the point
% n/2+1,n/2+1.

[n n1]=size(m);

p=ToPolar(m,n/2,n*3);
ps=mean(p,2);  % sum over thetas
ms=ToRect(ps,n);
