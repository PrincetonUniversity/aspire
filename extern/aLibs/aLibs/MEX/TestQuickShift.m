% TestQuickShift
% Check the speed of QuickShift compared to circshift

m=randn(30,30);
m2=m;
tic
for i=1:1000
    m=QuickShift(m,[3,42]);
end;
QuickShiftTime=toc

tic
for i=1:1000
    m2=circshift(m2,[3,42]);
end;
CircShiftTime=toc

diff=sum(sum(m2-m))

