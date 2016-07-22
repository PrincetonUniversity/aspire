function [bsp binctrs binwidths]=BinSpectrum(sp, width)
% given the 1d power spectrum, return a "binned" spectrum in which the
% higher-frequency points are averaged together.
if nargin<2
     width = 0.1;
end;
wconst=log(10);

n=numel(sp);
nbins=1;
binctrs(1)=0;   % proportional to frequency
binwidths(1)=1;
edge=2;
while edge<n
    nbins=nbins+1;
    w=floor(wconst*width*edge)+1;
    if edge+w>n
        w=n-edge+1;
    end;
    bsp(nbins)=mean(sp(edge:edge+w-1));
    binctrs(nbins)=round(edge-1+w/2);  % frequency point
    binwidths(nbins)=w;

    edge=edge+w;
end;
w=n-edge+1;
binctrs(nbins)=round(edge-1+w/2);
