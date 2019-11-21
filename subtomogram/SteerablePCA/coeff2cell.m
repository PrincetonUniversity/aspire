function [ Cell ] = coeff2cell( coeff )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
maxL = size(coeff,1);
Cell = cell(maxL,1);
for l = 1:maxL
    temp = [];
    for m = 1:2*l-1
        temp = cat(2,temp,coeff{l,(l-1)^2+m});
    end
    Cell{l} = temp;
end

end

