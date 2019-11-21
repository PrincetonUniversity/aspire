%{
INPUT:
    t: degree of representation
    alpha,beta,gamma: Euler angles
OUTPUT:
    Wt: degree t representation

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function Wt = realWignerD(t,alpha,beta,gamma)

Wt = wignerd(t,[alpha,beta,gamma]);

[T,Tinv] = realY_to_complexY(t);

Wt = real(Tinv*Wt*T);

end









