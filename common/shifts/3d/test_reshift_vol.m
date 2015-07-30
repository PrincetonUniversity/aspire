%
% Basic testing of reshift_vol.
%
% Yoel Shkolnisky, February 2014.

data=load('cleanrib');
vol=real(data.volref);

%% Odd-sized volume, intergral shift
vol1=vol;
vol1=vol1(1:end-1,1:end-1,1:end-1);
vol2=vol1(:,:,[5:end 1:4]);
dx=[0 0 4];
vol3=reshift_vol(vol2,-dx);  % First shift paramter is the x coordinate. 

figure(1);
clf;
view3d(vol1,1.0e-4,'b')
view3d(vol2,1.0e-4,'g')

figure(2);
clf;
view3d(vol1,1.0e-4,'b')
view3d(vol3,1.0e-4,'g')

err=norm(vol1(:)-vol3(:))/norm(vol1(:));
disp(err);


%% Odd-sized volume, arbitrary shift
vol1=vol;
vol1=vol1(1:end-1,1:end-1,1:end-1);
dx=[1.543,3.777, 3.123];
vol2=reshift_vol(vol1,dx);  
vol3=reshift_vol(vol2,-dx);  % First shift paramter is the x coordinate. 

figure(1);
clf;
view3d(vol1,1.0e-4,'b')
view3d(vol2,1.0e-4,'g')

figure(2);
clf;
view3d(vol1,1.0e-4,'b')
view3d(vol3,1.0e-4,'g')

err=norm(vol1(:)-vol3(:))/norm(vol1(:));
disp(err);

%% Even-sized volume, intergral shift
vol1=vol;
vol2=vol1(:,:,[5:end 1:4]);
dx=[0 0 4];
vol3=reshift_vol(vol2,-dx);  % First shift paramter is the x coordinate. 

figure(1);
clf;
view3d(vol1,1.0e-4,'b')
view3d(vol2,1.0e-4,'g')

figure(2);
clf;
view3d(vol1,1.0e-4,'b')
view3d(vol3,1.0e-4,'g')

err=norm(vol1(:)-vol3(:))/norm(vol1(:));
disp(err);

%% Even-sized volume, arbitrary shift
% Note that for even-sized volumes there are small errors when the shifts
% are not integral pixels (of the order of 10^-5).
vol1=vol;
dx=[1.543,3.777, 3.123];
vol2=reshift_vol(vol1,dx);  
vol3=reshift_vol(vol2,-dx);  % First shift paramter is the x coordinate. 

figure(1);
clf;
view3d(vol1,1.0e-4,'b')
view3d(vol2,1.0e-4,'g')

figure(2);
clf;
view3d(vol1,1.0e-4,'b')
view3d(vol3,1.0e-4,'g')

err=norm(vol1(:)-vol3(:))/norm(vol1(:));
disp(err);
