% function meq=ImEqualize(m,nbins)
% function meq=ImEqualize(m)
% Perform a histogram equalization of the image m.
% Return a transformed version of m where all values in the range 0...range
% are equally likely.  The transformation is monotonic.
% if nargin<2
%     nbins=2^16;
% end;
nbins=2^20;

sz=size(m);
nel=prod(sz);
range=256;
m=m-min(m(:)+eps); % make it positive.
mx=max(m(:));
mnorm=ceil(m(:)*(nbins-1)/mx); % values 1..nbins
sd=std(m(:));

h=hist(mnorm(:),1:nbins);
hc=cumsum(h);  % cumulative distribution

%%
qrange=(-range/2:range/2-1)*2/range;
kq=4;
qh=sign(qrange).*(exp(abs(kq*qrange))-1)/(exp(kq)-1);  % output range compression,
plot(qh)
%%
%  maps range -> range points with compression of extgremes
% meq=hc(max(1,mnorm))*range/nel;
mec=hc(max(1,mnorm))/nel; % uniform from 0..1
men=2*mec-1;  % -1 to 1
kq=2;
meq=sign(men).*(1-abs(men)).^kq; % -1 to 1
hist(meq,100);
meq=reshape((meq+1)*range/2,sz);
subplot(222); imac(meq);
subplot(221); hist(meq(:),100);