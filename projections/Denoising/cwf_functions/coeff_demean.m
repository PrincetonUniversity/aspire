function [ coeff ] = coeff_demean( data, R, basis, sample_points, num_pool )

% Compute Fb coefficients for images after mean subtraction
% data: Image stack (real space)
% basis: precomputed basis functions
% num_pool: Number of workers
% sample_points: Quadrature points
% Tejal: updated Oct 2015, Based on Jane's FFB code

n = size(data, 3);
data2 = cell(num_pool, 1);
nb = floor(n/num_pool);
remain = n - nb*num_pool;

for i = 1:remain
    data2{i} = data(:, :,(nb+1)*(i-1)+1: (nb+1)*i);
end;

count = (nb+1)*remain;
for i = remain+1:num_pool
    data2{i} = data(:, :, count + (i-remain-1)*nb+1: count + (i-remain)*nb);
end;

clear data;

[coeff ]= FBcoeff_nfft(data2, R, basis, sample_points, num_pool);


end

