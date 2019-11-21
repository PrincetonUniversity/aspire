function [ Cell ] = vec2cell( vec, siz_n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

maxL = size(siz_n,1)-1;
Cell = cell(maxL+1,1);
count = 0;
for ll = 0:maxL
    n = siz_n(ll+1);
    temp = [];
    for m = 1:2*ll+1
        count = count + n;
        temp = cat(2,temp,vec(count-n+1:count)); 
    end
    Cell{ll+1} = temp;
end

end

