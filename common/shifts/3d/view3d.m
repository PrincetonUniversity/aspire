function view3d(vol,threshold,clr)
%VIEW3D Render a 3D volume.
%   Render a 3D volume at the isosurface given by threshold.
%
%Yoel Shkolnisky, November 2013.

if ~exist('clr','var')
    clr='r';
end

p=patch(isosurface(vol,threshold));
isonormals(vol,p)
set(p, 'FaceColor',clr, 'EdgeColor','none')
daspect([1 1 1]);
view(3); 
axis vis3d tight;
box on; 
grid on;
camproj perspective

currentLight=findobj(gca,'Type','light');
if isempty(currentLight)
    camlight;
end
lighting phong, alpha(1)
axis([1 size(vol,1) 1 size(vol,2) 1 size(vol,3)])
xlabel('x'); ylabel('y'); zlabel('z');