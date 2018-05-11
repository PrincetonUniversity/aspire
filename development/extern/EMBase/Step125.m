function [xs unitval]=Step125(x)
% function [xs unitval]=Step125(x)
% Find the value in a 1..2..5 x 10^n sequence that is >= x.
% The returned unitval is the chosen value 1, 2 or 5.
%
n=floor(log10(x));
x0=x/10^n;
testvals=[0.5 1 2 5 10];
index=find((x0<=testvals),1);
unitval=testvals(index);
xs=unitval*10^n;
if unitval<1
    unitval=10*unitval;
end;