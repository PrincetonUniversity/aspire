function f=LeastCommonFactor(n1,n2)
% function f=LeastCommonFactor(n1,n2)
% find the product of all the common factors between two integers.
% e.g. LeastCommonFactor(2,16) returns 2.
%      LeastCommonFactor(768,1024) returns 256.
% For example to resize an image m from size n1 to n2, do this:
% lcf=LeastCommonFactor(n1,n2);
% imgOut=BinImage(ExpandImage(m,n2/lcf),n1/lcf);

f1=factor(n1);
f2=factor(n2);
f=1;
for i=1:numel(f1)
    fac=f1(i);
    q=find(f2==fac,1);
    if numel(q)>0
        f2(q)=1;
        f=f*fac;
    end;
end;
