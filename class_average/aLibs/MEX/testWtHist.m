% testWeightedHisto.m

n=1024;
r=Radius(n)+1;
% r=r(:);
r(1:100)=-1000;  % give some out-of-bounds points.
r(101:200)=1e6;
c=CTF(n,1,.02,4,2,50,0.1);
% c=Crop(c,n);
tic
[hist norm]=WeightedHisto(r,c,n);
toc
plot(hist./(norm+eps));
