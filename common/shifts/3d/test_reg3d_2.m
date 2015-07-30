% Plot the registration accuracy as a function of the noise level.
%
% Yoel Shkolnisky, January 2014.

%% Load the first volume
data=load('cleanrib');
vol1=real(data.volref);
vol1=vol1(1:end-1,1:end-1,1:end-1);

%% Create a shifted image
dx=[1.543,3.777, 3.123];
%dx=[1.01,2.02,3.03];
%dx=[1,2,3];
vol2=reshift_vol(vol1,dx);  % First shift paramter is the x coordinate. 

vol1=vol1./norm(vol1(:));
vol2=vol2./norm(vol2(:));

%% Create graph of registration accuray 
% For N different noise levels, compute the registration error for that 
% noise level, and show a graph of errors vs. sigma.

N=20;
sigma_arr=linspace(0,8,N);
err_arr=zeros(N,1);

Mfigs=ceil(sqrt(N)); % Number of rows in the subplot figure
Nfigs=ceil(N/Mfigs); % Number of columns in the subplot figure
h1=figure;

cc=floor(size(vol1,1)/2)+1; % Index of the center slice of the volume.

for j=1:N
    
    rng('default')
    rng(111);
    sigma=sigma_arr(j);
    n1=sigma*randn(size(vol1))./sqrt(numel(vol1(:)));
    n2=sigma*randn(size(vol2))./sqrt(numel(vol2(:)));
    vol1n=vol1+n1;
    vol2n=vol2+n2;

    % Show a slice of the noisy volume to assess the noise level.
    figure(h1);
    subplot(Mfigs,Nfigs,j);
    imagesc(vol1n(:,:,cc)); axis image; axis off; colormap(gray);
    title(sprintf('%5.3f',sigma));
    
    % Register the two volumes.
    [estdx,err]=register_translations_3d(vol1n,vol2n,dx);
    
    % Save the registration error.
    err_arr(j)=norm(err);
end

% Plot the results.
figure;
plot(sigma_arr,err_arr,'s-','LineWidth',2);
xlabel('sigma');
ylabel('Error (in pixels)');
grid on
