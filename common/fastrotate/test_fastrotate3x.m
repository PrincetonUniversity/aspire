function test_fastrotate3x
%TEST_FASTROTATE3X  Test the function fastrotate3x.
%
% The output of the function should be a small ball rotating around the
% x-axis CCW.
%
% Yoel Shkolnisky, November 2013.

vol=cryo_gaussian_phantom_3d('beadY_params',65,1);
figure;
camlight;
clrs=jet(36);
clridx=1;
for k=1:10:360
    volrx=fastrotate3x(vol,k);    
    view3d(volrx,0.5,clrs(clridx,:));
    clridx=clridx+1;
    pause(0.5)
    drawnow;
end
