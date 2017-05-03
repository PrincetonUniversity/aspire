% test_cryo_align_densities_C4_1
%
% Align two clean density maps which are C4 symmetric where one is
% rotate/translated relative to the other.
%
% Yoel Shkolnisky, March 2014

%initstate;

%% Generate C4 symmetric reference density map.
load cleanrib
vol=real(volref);
vol=GaussFilt(vol,0.8);
vol=(vol+fastrotate3z(vol,90)+fastrotate3z(vol,180)+fastrotate3z(vol,270))/4;

% Instead of the above lines you can use:
% vol=cryo_gaussian_phantom_3d('C1_params',64,1); % Perfect Gaussian, perfect results.

%% Generate a rotated copy of the reference.
[R,~,~]=svd(rand(3)); % Generate random rotation
%R=eye(3);
volRotated=fastrotate3d(vol,R); % Rotate the reference volume by the random rotation
volRotated=reshift_vol(volRotated,[5 0 0]);
%
% Note: it seems that the estimated shift estdx below is equation to (the
% negative of) R.'*[0 5 0].', where the flip in the position of the 5 is
% due to the difference convetions of which coordinate is the x axis
% between the difference functions. Check if the estdx should R.'*[0 5 0].'
% or R*[0 5 0].'

%% Visualize the two maps.
figure(1); clf; view3d(vol); title('Reference volume');
figure(2); clf; view3d(volRotated); title('Rotated volume');

%% Align
verbose=1;
tic;
[Rest,estdx,vol2aligned]=cryo_align_densities_C4(vol,volRotated,0,verbose,R);
toc

figure(3); clf; view3d(vol2aligned,2.0e-4); title('Aligned volume');

% Correlation between two original volumes
c1=corr(vol(:),volRotated(:));
fprintf('Correlation between two original volumes %7.4f\n',c1);

% Correlation after alignment
c2=corr(vol(:),vol2aligned(:));
fprintf('Correlation between original and aligned volume %7.4f\n',c2);

% The original and aligned volume should have high agreement according to
% the FSC. FSC may degrade due to sub-pixel misalignment errors.
plotFSC(vol,vol2aligned,0.5,1);