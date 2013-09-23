% Properly sampling a Gaussian in two dimensions with an even number of
% points. The points are distributed symmetrically around 0.
%
% Yoel Shkolnisky, September 2013.

sig2 = (0.1)^2;  % Sigma squared.

g = @(x) (1/(2*pi*sig2))*exp(-sum(x.^2,2)/(2*sig2));  % A Gaussian
fg = @(omega) exp(-sum(omega.^2,2)*sig2/2); % It's Fourier transform

T=20*sqrt(sig2);  % Sampling interval for the Gaussian
N=128;   % Number of points.
dx=T/N;  % Spacing between sampling points.
dw=2*pi/T;  % Spacing between FFT sampling points (in frequnecy domain)
K=-N/2:N/2-1; 
J=-N/2:N/2-1;

xk=(K+1/2)*dx; % Sampling points in spatial domain
[Xk,Yk]=ndgrid(xk,xk);
wj=(J+1/2)*dw; % Sampling points in frequnecy domain
[Wx, Wy]=ndgrid(wj,wj);

gk=g([Xk(:) Yk(:)]);      % Evaluate the Gaussian
gk=reshape(gk,N,N);
imagesc(gk);

hatg = cfft2e(gk); % Take centered FFT
disp(norm(imag(hatg))/norm(hatg)) % Should be of the order of machine precision
hatg=real(hatg); % Get rid of tiny imaginary components
hatg=hatg.*(dx)^2; % Normalize to get a Fourier transform

% Compare to analytical Fourier transform
hatgref=fg([Wx(:) Wy(:)]);
hatgref=reshape(hatgref,N,N);
disp(norm(hatgref-hatg)/norm(hatgref)) % Should be of the order of machine precision.

% Now test icfft2e.
% Make sure that icfft2e recovers the original image.
hatg=hatg./(dx.^2); % Get rid of the normalization;
gk2=icfft2e(hatg);
disp(norm(gk(:)-gk2(:))/norm(gk(:))); % Should be of the order of machine precision.
