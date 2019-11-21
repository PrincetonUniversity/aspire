function [ F ] = coeff2real( coeff, basis, sample_points )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ dataS ]= FBcoeff_sphere_inv(coeff, basis, sample_points);
fS = [];
for i = 1:n_r
    tempS = ssht_inverse(dataS(i,:),L);
    fS = cat(1,fS,tempS(:));
end



end

