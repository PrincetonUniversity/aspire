% Similar to test_CenterOfMass_3 but project in a random direction and not
% in the z direction.
%
% Tests both even and odd sized volume/image.
%
% Yoel Shkolnisky, November 2014.

%% Odd sized test
initstate;

clear;
voldef='C1_params';  % A simple centered phantom.
n=65;               % Size of the volume.
rmax=1;              % The volume is generated by sampling the phantom on 
                     % a 3D Cartesian grid where each dimension has n
                     % samples on [-rmax,rmax]
vol=cryo_gaussian_phantom_3d(voldef,n,rmax);  % Generate the nxnxn volume.
q=qrand(1);

% First center the volume then project.
cmv=CenterOfMass(vol);
vol_aligned=reshift_vol(vol,cmv);
proj1=cryo_project(vol_aligned,q,n,'double'); % Double precision accuracy 
                                              % so we don't introduce more
                                              % errors 

% First project then center
proj=cryo_project(vol,q,n,'double');
cmp=CenterOfMass(proj);
proj2=reshift_image(proj,cmp);

% Due to the shape of the phantom and the projection in an arbitrary
% direction, the error is not very small. For a simpler phantom such as
% one_ball the error will be double precision.
err=norm(proj1(:)-proj2(:))/norm(proj1(:));
fprintf('Error odd %6.4e\n',err);
assert(err<1.0e-4);

%% Even sized test
clear;
voldef='C1_params';  % A simple centered phantom.
n=64;               % Size of the volume.
rmax=1;              % The volume is generated by sampling the phantom on 
                     % a 3D Cartesian grid where each dimension has n
                     % samples on [-rmax,rmax]
vol=cryo_gaussian_phantom_3d(voldef,n,rmax);  % Generate the nxnxn volume.
q=qrand(1);

% First center the volume then project.
cmv=CenterOfMass(vol);
vol_aligned=reshift_vol(vol,cmv);
proj1=cryo_project(vol_aligned,q,n,'double'); % Double precision accuracy 
                                              % so we don't introduce more
                                              % errors 

% First project then center
proj=cryo_project(vol,q,n,'double');
cmp=CenterOfMass(proj);
proj2=reshift_image(proj,cmp);

% Due to the shape of the phantom and the projection in an arbitrary
% direction, the error is not very small. For a simpler phantom such as
% one_ball the error will be double precision.
err=norm(proj1(:)-proj2(:))/norm(proj1(:));
fprintf('Error even %6.4e\n',err);
assert(err<1.0e-4);

%% Print OK.
fprintf('Test OK\n'); % No assertions have been violated.
