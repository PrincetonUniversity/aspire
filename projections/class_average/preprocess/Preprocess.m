function [ data ] = Preprocess (data, defocus_group, c)
%This function estimate noise power spectrum from the non-particle region
%   Input: data
%          defocus_group: the defocus group indices for each projection
%          c: estimated 2D CTFs
%   Output: data
%   Zhizhen Zhao June 2013

n=size(data, 3);
mean_data=mean(data, 3);
for i=1:n
    data(:, :, i)=data(:, :, i)-mean_data;
end;
[ P ] = Noise_Estimation( data );
[ P2 ] = radial_average( P );
[ data ] = cryo_prewhiten(data, P2);
[ data ] = Phase_Flip( data, defocus_group, c );

end
