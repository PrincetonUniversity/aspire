function [ eul ] = ind2euler( ind, B )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
k = floor(ind/(2*B)^2);
j1 = floor(mod(ind,(2*B)^2)/(2*B));
j2 = mod(mod(ind,(2*B)^2),2*B)-1;
alpha = 180*j1/B;
beta = 180*(2*k+1)/(4*B);
gamma = 180*j2/B;
eul = [alpha beta gamma];
end

