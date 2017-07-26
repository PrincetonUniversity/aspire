% CRYO_PERIODOGRAM_MEAN Calculate mean of periodograms
%
% Usage
%    x_per_mean = cryo_periodogram_mean(x_per);
%
% Input
%    x_per: An array of periodograms of size L1-by-L2-by-n obtained from
%       cryo_periodogram.
%
% Output
%    x_per_mean: The mean of the periodograms.

function x_per_mean = cryo_periodogram_mean(x_per)
    n = size(x_per, 3);

    x_per_mean = sum(x_per, 3)/n;
end
