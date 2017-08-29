% Test the function cryo_add shifts.
%
% Yoel Shkolnisky, September 2013.

%% Test 1: Shift back and forth
% Shift a projection and then shift it back.
clear;
voldef='C1_params';  % Can be 'one_ball' for a simpler phantom.
rot = eye(3);        % Identity rotation - projection along the z direction.
n=65;                % Size of the projections. Can be even or odd.
                     % Sampling must be sufficient for the projected
                     % (Guassian) phantom.
rmax=1;              % The value of n must be large enough for the given 
                     % rmax. See "Test 1" for more details.
p1=cryo_project_gaussian(voldef,n,rmax,rot); % Compute analytic projection.
sp1=cryo_addshifts(p1,[20 0]);   % Shift 
sp2=cryo_addshifts(sp1,[-20 0]); % and Shift back
err=norm(p1(:)-sp2(:))/norm(p1(:));
disp(err); % Error should be close to machine precision.
subplot(1,3,1);imagesc(p1); colorbar; axis image;  % Original
subplot(1,3,2);imagesc(sp2); colorbar; axis image; % Shifted
subplot(1,3,3); imagesc(p1-sp2); colorbar; axis image; % Difference

%% Test 2: Sub-pixel shift
% Shift a projection by a non-integer shift. Inspect the results visually.
clear;
voldef='C1_params';  % Can be 'one_ball' for a simpler phantom.
rot = eye(3);        % Identity rotation - projection along the z direction.
n=65;                % Size of the projections. Can be even or odd.
                     % Sampling must be sufficient for the projected
                     % (Guassian) phantom.
rmax=1;              % The value of n must be large enough for the given 
                     % rmax. See "Test 1" for more details.
p1=cryo_project_gaussian(voldef,n,rmax,rot); % Compute analytic projection.
sp1=cryo_addshifts(p1,[10.5 0]);   % Shift 
subplot(1,2,1);imagesc(p1); colorbar; axis image;  % Original
subplot(1,2,2);imagesc(sp1); colorbar; axis image; % Shifted

%% Test 3: Random shift
% Shift the projection by a random shift.
clear;
voldef='C1_params';  % Can be 'one_ball' for a simpler phantom.
initstate;
rot = rand_rots(1); % Identity rotation - projection along the z direction.
n=65;                % Size of the projections. Can be even or odd.
                     % Sampling must be sufficient for the projected
                     % (Guassian) phantom.
rmax=1;              % The value of n must be large enough for the given 
                     % rmax. See "Test 1" for more details.
p1=cryo_project_gaussian(voldef,n,rmax,rot); % Compute analytic projection.
[sp1,ref_shifts]=cryo_addshifts(p1,[],10,5);   % Shift 
subplot(1,2,1);imagesc(p1); colorbar; axis image;  % Original
subplot(1,2,2);imagesc(sp1); colorbar; axis image; % Shifted
fprintf('Shifts are [%d %d]\n',ref_shifts(1),ref_shifts(2));

%% Test 4: Compare to gen_projections_v2
initstate;
[P1,~,ref_shifts,rots]=cryo_gen_projections(1,1,5,1);
n=size(P1,1);
%volref=cryo_gaussian_phantom_3d('C1_params',n,1);
load cleanrib
P2=cryo_project(volref,rots,n);
P2=permute(P2,[2 1]);
SP2=cryo_addshifts(P2,ref_shifts);
subplot(1,3,1); imagesc(P1); colorbar; axis image;
subplot(1,3,2); imagesc(SP2); colorbar; axis image;
subplot(1,3,3); imagesc(P1-SP2); colorbar; axis image;
err=norm(P1(:)-SP2(:))/norm(P1(:));
disp(err); % Should be of the order of machine precision.

