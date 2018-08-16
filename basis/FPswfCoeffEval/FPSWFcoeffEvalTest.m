%clc;
clear variables;
close all;

%% Define parameters
beta = 1;
L = 128;
c = beta*pi*L;
T = 1e1;
realFlag = 0;
numChunks = 64;

% - Choose tolerances
phiApproxErr = 1e-16;
radialQuadErr = 1e-16;
angularQuadErr = 1e-16;

%% Generate quadrature nodes and weights
% [quadRulePtsX,quadRulePtsY,quadRuleWts,radialQuadPts,quadRuleRadialWts,numAngularPts] = generatePswfQuad(4*L,2*c,phiApproxErr,radialQuadErr,angularQuadErr,realFlag);

%% Generate PSWFs on the quadrature nodes
% [ PSWF_mat_quad, Alpha_Nn, ~, ~ ] = PSWF_2D_gen_with_truncation( L, beta, eps, T, L*quadRulePtsX.', L*quadRulePtsY.' );

%% Test orthonormality of PSWFs on disk using quadratures
% orthonormalityTest = PSWF_mat_quad' * diag(quadRuleWts) * PSWF_mat_quad;
% max(max(abs(triu(orthonormalityTest,1))))

%% Generate rectangular grids 
% x_1d_grid = -1:(1/L):1;   % - Odd number of points
x_1d_grid = -1:(1/L):(1-1/L);   % - Even number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= 1);

x_cartesian = x_2d_grid(points_inside_the_circle);
y_cartesian = y_2d_grid(points_inside_the_circle);

%% Prepare quadratures-based transform and evaluate spectrum 
% USFFT_mat = exp(-1i*c*[quadRulePtsX.' quadRulePtsY.'] * [x_cartesian.'; y_cartesian.']);
% G = c/(2*pi*L) * PSWF_mat_quad' * bsxfun(@times,quadRuleWts.',USFFT_mat);
% 
% spectrumTest = G*G';
% s = eig(spectrumTest);
% figure; bar(abs(s));

%% generate images
nImages = 2e4;
images = randn(2*L,2*L,nImages);

%% Compute coefficients using standard inner product and also using quadratures
[ coeffVecQuadFast, normalizedPSWF_mat] = fastPswfCoeffEval( images, x_cartesian, y_cartesian, L, beta, T, numChunks, realFlag );

images = reshape(images,size(images,1)*size(images,2),size(images,3));
images = images(points_inside_the_circle(:),:);
coeffVecStandard = normalizedPSWF_mat' * images;

% coeffVecQuad = G * images;

%% Compare resulting expansion coefficients of different methods
% err1 = max(max(abs(coeffVecQuadFast - coeffVecQuad)))
% err2 = max(max(abs(coeffVecQuadFast - coeffVecStandard)))
err2 = mean(max(abs(coeffVecQuadFast - coeffVecStandard)))

