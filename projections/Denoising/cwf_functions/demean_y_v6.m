function [data] = demean_y_v6(data, F, EstMean, index)

% Subtract estimated "mean image" from stack of images
% y - A \mu : Demeaned y
% New y_i = y_i - A_i \mu
% Do in place
% Tejal Bhamre, Oct 2015

for i=1:size(data,3)
    data(:,:,i) = data(:,:,i)-(F(:,:,index(i)).*EstMean);
end


end
