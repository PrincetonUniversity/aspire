function test_fastrotate3z
%TEST_FASTROTATE3Z  Test the function fastrotate3z.
%
% The output of the function should be a small ball rotating around the
% z-axis CCW.
%
% Yoel Shkolnisky, November 2013.

vol=cryo_gaussian_phantom_3d('beadX_params',65,1);
figure;
camlight;
clrs=jet(36);
clridx=1;
for k=1:10:360
    volrx=fastrotate3z(vol,k);
    view3d(volrx,0.5,clrs(clridx,:));
    clridx=clridx+1;
    pause(0.5)
    drawnow;
end
