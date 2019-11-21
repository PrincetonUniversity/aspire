%{
INPUT:
    t: degree of representation
OUTPUT:
    T: real Ylm to complex Ylm
    Tinv: complex Ylm to real Ylm

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ T , Tinv ] = realY_to_complexY_noJ(t)

Tinv = zeros(2*t+1);
for k = (-t):t
    
    if k < 0
        Tinv(k+t+1,k+t+1) = 1i/sqrt(2);
        Tinv(k+t+1,-k+t+1) = -1i*(-1)^k/sqrt(2);
    elseif k == 0
        Tinv(t+1,t+1) = 1;
    else
        Tinv(k+t+1,k+t+1) = (-1)^k/sqrt(2);
        Tinv(k+t+1,-k+t+1) = 1/sqrt(2);
    end
    
end

T = zeros(2*t+1);
for k = (-t):t
    
    if k < 0
        T(k+t+1,k+t+1) = -1i/sqrt(2);
        T(k+t+1,-k+t+1) = 1/sqrt(2);
    elseif k == 0
        T(k+t+1,k+t+1) = 1;
    else
        T(k+t+1,k+t+1) = (-1)^k/sqrt(2);
        T(k+t+1,-k+t+1) = 1i*(-1)^k/sqrt(2);
    end
    
end

end









