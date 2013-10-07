function out=MirrorX(in)
% Flip a 1D, 2D or 3D image about the X-axis
% Works correctly for FFT origin.

nx=size(in,1);
if mod(nx,2)>0  % it's odd, so it's simple
    newX=nx:-1:1;
else
    newX=[1 nx:-1:2]; % leave the x=1 column alone (periodic bdry)
end;
out=in(newX,:,:);
