function [hist,norm]=WtHist(indices, weights, nbins)
% function [hist,norm]=WtHist(indices, weights, nbins)
% M-file implementation of the weighted-histogram MEX file.
% The MEX function runs about 3x faster.

hist=zeros(nbins,1);
norm=zeros(nbins,1);
if nbins < 1
    return
end;

for i=1:numel(indices)
    k=min(nbins,max(1,indices(i)));
    hist(k)=hist(k)+weights(i);
    norm(k)=norm(k)+1;
end;
