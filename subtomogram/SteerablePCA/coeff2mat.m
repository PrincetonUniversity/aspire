function [ A ] = coeff2mat( coeff )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

l = size(coeff,1);
A = zeros(length(coeff{1,1}),l^2);
for i = 1:l
    for j = 1:2*i-1
        len = length(coeff{i,(i-1)^2+j});
        A(1:len,(i-1)^2+j) = coeff{i,(i-1)^2+j};
    end
end


end

