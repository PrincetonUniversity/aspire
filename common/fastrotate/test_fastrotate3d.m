%TEST_FASTROTATE3D Test the function fastrotate3d.
%   
%Yoel Shkolnisky, November 2013.

%% First step - measure error for each of the three basic rotations.
% This shows that the rotation functions are accurate.

%% Measure error of rotate3x
vol=cryo_gaussian_phantom_3d('beadY_params',65,1); % Use N=129 to get an error of 10^-12.
volr=vol;
% Rotate by 1 degree 360 times and check error.
figure;
for k=1:360
    volr=fastrotate3x(volr,1);
    if mod(k,10)==0
        view3d(volr,2.05e-4);        
        drawnow;
    end
end
clf;
subplot(1,2,1);
view3d(vol,2.05e-4);
subplot(1,2,2);
view3d(volr,2.05e-4);

% Display error
err=norm(vol(:)-volr(:))/norm(vol(:));
fprintf('Error in rotation: %e\n',err);

%% Measure error of rotate3y
vol=cryo_gaussian_phantom_3d('beadX_params',65,1); % Use N=129 to get an error of 10^-12.
volr=vol;
% Rotate by 1 degree 360 times and check error.
figure;
for k=1:360
    volr=fastrotate3y(volr,1);
    if mod(k,10)==0
        view3d(volr,2.05e-4);        
        drawnow;
    end
end
subplot(1,2,1);
view3d(vol,2.05e-4);
subplot(1,2,2);
view3d(volr,2.05e-4);

% Display error
err=norm(vol(:)-volr(:))/norm(vol(:));
fprintf('Error in rotation: %e\n',err);

%% Measure error of rotate3z
vol=cryo_gaussian_phantom_3d('beadZ_params',65,1); % Use N=129 to get an error of 10^-12.
volr=vol;
% Rotate by 1 degree 360 times and check error.
figure;
for k=1:360
    volr=fastrotate3z(volr,1);
    if mod(k,10)==0
        view3d(volr,2.05e-4);        
        drawnow;
    end
end
subplot(1,2,1);
view3d(vol,2.05e-4);
subplot(1,2,2);
view3d(volr,2.05e-4);

% Display error
err=norm(vol(:)-volr(:))/norm(vol(:));
fprintf('Error in rotation: %e\n',err);

%% Check that the coordinates axes are rotated properly.
% Applying the rotation R on a volume containing a bead at location [1 0
% 0].' is equal to applying R on the vector [1 0 0].'. Since the latter is
% equal to the first column of R, the bead should be transformed to the
% point in R^{3} whose coordinates are the first column of R. The same
% holds for the second column of R and a bead at [0 1 0].', and for the
% third column of R and a bead at [0 0 1].'.
%
% Note that the vector R1, R2, and R3 are not precisely equal to the
% corresponding columns of the matrix R, since the point to which a bead
% has been rotated is found using [vv,ii]=max(...), whose accuracy is not
% high since is likely that the center of a bead won't coincide with one of
% the grid points. Since the Gaussian of the bead decays fast, deviations
% can be rather large. Thus difference between the two matrices computed
% below is abour 10^-2. However, the three column of the two matrices are
% clearly in the same direction, which means that the function fastrotate3d
% applies the matrix R to the input volume.

initstate;
N=65;
R=rand_rots(1);

volx=cryo_gaussian_phantom_3d('beadX_params',N,1);
volrx=fastrotate3d(volx,R);
[vv,ii]=max(volrx(:));  % Find the center of the rotated volume.
[y,x,z]=ind2sub([N N N],ii);  % Not that y and x are switched, since in matlab the first coordinate of an image is actually the y axis.
% The bead is at [1/2 0 0]. The following line transforms the center of the
% detected volume into [-1/2 1/2], then multiplies by 2 to strech it [-1 1]
% and then multiplies by another two since the original volume was centeted
% at [1/2 0 0] and we want to transform it to [1 0 0]
R1=([x y z]-(N+1)/2)/N*4; 

% Same for the y axis
voly=cryo_gaussian_phantom_3d('beadY_params',N,1);
volry=fastrotate3d(voly,R);
[vv,ii]=max(volry(:));  
[y,x,z]=ind2sub([N N N],ii);  
R2=([x y z]-(N+1)/2)/N*4; 

% Same for the z axis
volz=cryo_gaussian_phantom_3d('beadZ_params',N,1);
volrz=fastrotate3d(volz,R);
[vv,ii]=max(volrz(:));  
[y,x,z]=ind2sub([N N N],ii);  
R3=([x y z]-(N+1)/2)/N*4; 
Rest=[R1 ; R2 ; R3].';

disp('Original:');
disp(R)
disp('Estimated from rotated volume:');
disp(Rest);
fprintf('err=%e\n',norm(R-Rest,'fro'));


%% Finally, test that rotating th volume back and forth restres the original volume.
vol=cryo_gaussian_phantom_3d('C1_params',65,1); % Use N=129 for double precision accuracy.
R=rand_rots(1);

[psi,theta,phi]=rot2xyz(R.');
assert(norm((Rz(phi)*Ry(theta)*Rx(psi))*R-eye(3))<1.0e-13)
volr=fastrotate3d(vol,R);
volir=fastrotate3d(volr,R.');

figure;
subplot(1,3,1);
view3d(vol,2.05e-4);
subplot(1,3,2);
view3d(volr,2.05e-4);
subplot(1,3,3);
view3d(volir,2.05e-4);
norm(vol(:)-volir(:))/norm(vol(:))
