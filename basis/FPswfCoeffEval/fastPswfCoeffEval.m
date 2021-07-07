function [ coeffVecQuadFast, normalizedPSWF_mat, Alpha_Nn, ang_freqs, rad_freqs, FPSWF_coeff_time ] = fastPswfCoeffEval( images, x, y, L, beta, T, numChunks, realFlag )
% 
% addpath(genpath('/opt/nfft-3.2.3'));
% addpath(genpath('/home/yoel/Downloads/nfft-3.3.2'));

nImages = size(images,3);
images = reshape(images,size(images,1)*size(images,2),nImages); % Reshape from 3d array to a matrix

c = beta*pi*L;

% - Choose tolerances
phiApproxErr = 1e-16;
radialQuadErr = 1e-16;
angularQuadErr = 1e-16;

%% Generate quadrature nodes and weights
[quadRulePtsX,quadRulePtsY,~,radialQuadPts,quadRuleRadialWts,numAngularPts] = generatePswfQuad(4*L,2*c,phiApproxErr,radialQuadErr,angularQuadErr,realFlag);

%% Generate rectangular grids 
x_1d_grid = -1:(1/L):1;   % - Odd number of points
%% x_1d_grid = -1:(1/L):(1-1/L);   % - Even number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= 1);
image_height = numel(x_1d_grid);
points_inside_the_circle_vec = reshape(points_inside_the_circle, image_height*image_height, 1);

%% Take only samples inside the unite disk
%% images = images(points_inside_the_circle(:),:);
images = images(points_inside_the_circle_vec(:),:);

%% Generate PSWF basis on Cartesian grid 
[ normalizedPSWF_mat, Alpha_Nn, ang_freqs, rad_freqs ] = PSWF_basis_gen_v3( L+1, L, beta, phiApproxErr, T, L*x, L*y );

% normalizedPSWF_mat = 0;
% Alpha_Nn = 0;
% ang_freqs = 0;
% rad_freqs = 0;

%% Evaluate radial part of PSWFs on quadrature nodes
[ PSWF_radial_quad, ~, ang_freq, ~] = PSWF_2D_radial_quad( L, beta, phiApproxErr, T, radialQuadPts);

%% Compute coefficients using fast quadrature-based method
tic
% - Partition the images into 'numChunks' chunks, and apply NFFT to them in parallel
chunkSize = ceil(nImages/numChunks);
usFftPts = c/L*[quadRulePtsX.' quadRulePtsY.'];
% - prepare chunks
imagesCell = cell(1,numChunks);
coeffVecQuadFast_cell = cell(1,numChunks);
for i = 1:numChunks
    imagesCell{i} = images(:,(chunkSize*(i-1)+1):min(chunkSize*i,nImages));
end
% - Compute NFFT and expansion coefficients
parfor i = 1:numChunks
    nfftRes = computeNfft( imagesCell{i}, usFftPts, L, points_inside_the_circle_vec);
    coeffVecQuadFast_cell{i} = fastPswfIntegration(nfftRes ,c ,L, numAngularPts, ang_freq, radialQuadPts, quadRuleRadialWts, PSWF_radial_quad, realFlag);
    imagesCell{i} = [];
end

% - Gather coefficients
coeffVecQuadFast = zeros(numel(ang_freq),nImages);
for i = 1:numChunks
    coeffVecQuadFast(:,(chunkSize*(i-1)+1):min(chunkSize*i,nImages)) = coeffVecQuadFast_cell{i};
    coeffVecQuadFast_cell{i} = [];
end

FPSWF_coeff_time = toc;
% display(['Computing the expansion coefficients took ',num2str(FPSWF_coeff_time),' seconds']);

end

