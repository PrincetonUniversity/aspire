function [ coeff_pos_k ] = cell2coeff( Cell )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
maxL = size(Cell,1);
coeff_pos_k = cell(maxL,maxL^2);
for l = 1:maxL
    coeff = Cell{l};
    for m = 1:2*l-1
        coeff_pos_k{l,(l-1)^2+m} = coeff(:,m);
    end
end

end

