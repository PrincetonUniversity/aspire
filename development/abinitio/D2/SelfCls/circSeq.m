%Make a circular sequence of integers between n1 and n2 modulo L
%Precond: abs(n1-n2)<=L
function [seq]=circSeq(n1,n2,L)
if abs(n1-n2)>L
    error('abs(n1-n2)>L');
end
if n2<n1
    n2=n2+L;
end
seq=mod(n1:n2,L);
seq(abs(seq)<1e-08)=L;

