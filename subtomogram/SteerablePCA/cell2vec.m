function [ vec ] = cell2vec( Cell )
% convert cell structured SB coefficient to a column vector, where n
% changes fastest, then m, then l
%   Detailed explanation goes here

maxL = length(Cell);
vec = [];
for ll = 1:maxL
    vec = cat(1,vec, Cell{ll}(:));
end

end

