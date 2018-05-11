function m1=RemoveGradients(m)
% function m1=RemoveGradients(m)
% As in CTFFIND3, remove linear gradients in 2D from an image.

[n1 n2]=size(m);
[X Y]=ndgrid(0:1/(n1-1):1,0:1/(n2-1):1);

mx0=mean(m(1,:));
mx1=mean(m(n1,:));
my0=mean(m(:,1));
my1=mean(m(:,n2));

m1=m-mx0-my0-(mx1-mx0)*X-(my1-my0)*Y;
