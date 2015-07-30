% Test the function register_translations_3d
%
% Yoel Shkolnisky, January 2014.

%% Load a volume
data=load('cleanrib');
vol1=real(data.volref);

% For odd sized volumes the answer is much more accurate. This can probably
% because in the even case there is a frequnecy without its conjugate
% conuterpart. This may be sloved by using frequencies that are centered
% around zero and using cffte and icffte.
vol1=vol1(1:end-1,1:end-1,1:end-1); 
%% Create a shifted volume
dx=[1.543,3.777, 3.123];
%dx=[1.001,2.002,3.003];
%dx=[1,2,3];
vol2=reshift_vol(vol1,dx);  % First shift paramter is the x coordinate. 

% Show the volme and its shifted copy
view3d(vol1,1.0e-4,'b')
view3d(vol2,1.0e-4,'g')

%% Add noise to the images.
sigma=0.1; % If sigma=0 (no noise), then the shift dx should be recovered very accurately.
vol1=vol1./norm(vol1(:));
vol2=vol2./norm(vol2(:));

rng('default')
rng(111);
n1=sigma*randn(size(vol1))./sqrt(numel(vol1(:)));
n2=sigma*randn(size(vol2))./sqrt(numel(vol2(:)));
vol1=vol1+n1;
vol2=vol2+n2;

%% Register
estdx=register_translations_3d(vol1,vol2,dx);

% Check error
vol3=reshift_vol(vol2,estdx);
disp(norm(vol3(:)-vol1(:))/norm(vol1(:)));
