function [ powerSpec , c , R ] = getDataParameters( varNoise )

n0 = 10000;
max_shift = 0;
step_size = 1;
noise_type = 'gaussian';
silent = 1;

[ projections , ~ , ~ , ~ ] = gen_projections( n0 , 1000 , max_shift , step_size , noise_type , silent );
if size(projections,1) ~= size(projections,2)
    error('projections not square')
end
L = size(projections,1);

noisy_projections = zeros(size(projections));
for i = 1:size(projections,3)
    noisy_projections(:,:,i) = projections(:,:,i) + sqrt(varNoise)*randn(L);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate bandlimit and compact support radius

[ c , R ] = avg_pspec( noisy_projections , varNoise );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate power spectrum

% remove mean from the img
mean_data = mean(noisy_projections,3);
noisy_projections = bsxfun( @minus , noisy_projections , mean_data );

% compute the power spectrum
for i = 1:size(noisy_projections,3)
    noisy_projections(:,:,i) = abs( cfft2( noisy_projections(:,:,i) ) ).^2;
end
powerSpec = mean(noisy_projections,3);

end









