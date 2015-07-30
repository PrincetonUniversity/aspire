% Plot the registration accuracy as a function of the noise level.
%
% Yoel Shkolnisky, January 2014.

%% Load the first image
data=load('projs');
im1=data.projs(:,:,10);

%% Create a shifted image
dx=[1.543,3.777]; % The shift to recover.
%dx=[2,0];
im2=reshift_image(im1,dx);  % First shift paramter is the x coordinate.
im1=im1./norm(im1(:));
im2=im2./norm(im2(:));


%% Create graph of registration accuray 
% For N different noise levels, compute the registration error for that 
% noise level, and show a graph of errors vs. sigma.

N=20;
sigma_arr=linspace(0,4,N);
err_arr=zeros(N,1);

Mfigs=ceil(sqrt(N)); % Number of rows in the subplot figure
Nfigs=ceil(N/Mfigs); % Number of columns in the subplot figure
h1=figure;

for j=1:N
    
    % Add noise to im1 and im2
    rng('default')
    rng(834);
    sigma=sigma_arr(j);
    n1=sigma*randn(size(im1))./sqrt(numel(im1(:)));
    n2=sigma*randn(size(im2))./sqrt(numel(im2(:)));
    im1n=im1+n1;
    im2n=im2+n2;

    % Show one of the noisy images to assess the noise level.
    figure(h1);
    subplot(Mfigs,Nfigs,j);
    imagesc(im1n); axis image; axis off; colormap(gray);
    title(sprintf('%5.3f',sigma));
    
    % Register the two images.
    [estdx,err]=register_translations_2d(im1n,im2n,dx);
    
    % Save the registration error
    err_arr(j)=norm(err);
end

% Plot the results.
figure;
plot(sigma_arr,err_arr,'s-','LineWidth',2);
xlabel('sigma');
ylabel('Error (in pixels)');
grid on
