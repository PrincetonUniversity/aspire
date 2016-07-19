function [ c_limit, R_limit ] = avg_pspec(img, var_hat, energy_thresh)
%Description
%Estimate bandlimit and compact support radius from images
%	Input:
%		img: image stack of size L x L x n
%		var_hat: estimated noise variance
%	Output:
%		c_limit: bandlimit.
%		R_limit: essential compact support radius.
%

L = size(img, 1);
if mod(L, 2) == 1
   N = floor(L/2);
   [ x, y ] = meshgrid( -N:N, -N:N);
else
   N = L/2;
   [ x, y ] = meshgrid( -N: N-1, -N:N-1);
end;
r = sqrt(x.^2 + y.^2);
mean_data = mean(img, 3);
%remove mean from the data
img = bsxfun(@minus, img, mean_data);

variance_map = var(img, [], 3);
%mean 2D variance radial function
radial_var = zeros(N, 1);
for i = 1:N
    radial_var(i) = mean(variance_map(r>=i-1 & r<i));
end;

%for i = 1:size(img, 3)
%    img(:, :, i) = abs(cfft2(img(:, :, i))).^2;
%end;  

%Tejal 
img=abs(cfft2(img)).^2;

pspec = mean(img, 3);
radial_pspec = zeros(N, 1);
%compute the radial power spectrum;
for i = 1:N
    radial_pspec(i) = mean(pspec(r>=i-1 & r<i))/L^2;
end;
%subtract the noise variance
%var_hat = 1;
%save tmp radial_pspec radial_var
radial_pspec = radial_pspec-var_hat;
radial_var = radial_var-var_hat;

% compute the cumulative variance and power spectrum.
c = 0.5/N:0.5/N:0.5;
c = c';
R = 0:N-1;
R = R';
cum_pspec = zeros(N, 1);
cum_var = zeros(N, 1);
for i = 1:N
    cum_pspec(i) = sum(radial_pspec(1:i).*c(1:i));
    cum_var(i) = sum(radial_var(1:i).*R(1:i));
end;
cum_pspec = cum_pspec/cum_pspec(end);
cum_var = cum_var/cum_var(end);

plot(cum_pspec);
figure; plot(cum_var);

c_limit = c(find(cum_pspec>energy_thresh, 1)-1);
R_limit = R(find(cum_var>energy_thresh, 1) - 1);

end
