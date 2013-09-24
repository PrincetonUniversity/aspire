function fwt = writeMRC(vol, pixelsize, filename);
% function fwt = writeMRC(vol, pixelsize, filename): write a 3D map into the mrc
% format (little-endian assumed).
fwt=0;
% Here is where we specify little-endian data:
fd=fopen(filename,'wb','ieee-le');
[nx, ny, nz]=size(vol);

% the first 10 values, which are integers:
% nc nr ns mode ncstart nrstart nsstart nx ny nz
% MODE     data type :
%      0        image : signed 8-bit bytes range -128 to 127
%      1        image : 16-bit halfwords
%      2        image : 32-bit reals
%      3        transform : complex 16-bit integers
%      4        transform : complex 32-bit reals
mode=2;
a(1:10)=[nx,ny,nz,mode,nx/2+1,ny/2+1,nz/2+1, nx-1,ny-1, nz-1];
cnt=fwrite(fd,a,'int32');
% a(1:10)

% the next 12, which are floats.
% the first three are the cell dimensions.
% xlength ylength zlength alpha beta gamma mapc mapr maps
% amin amax amean.
mapc=1;
mapr=2;
maps=3;
amin=min(min(min(vol)));
amax=max(max(max(vol)));
amean=mean(mean(mean(vol)));
b=[nx*pixelsize,ny*pixelsize,nz*pixelsize,ones(1,3)*90,1,2,3,amin,amax,amean];
cnt=fwrite(fd,b,'float32');
% b
%rez=b(1)/a(1);
% rez

%the next 30, which brings us up to entry no. 52.
cnt=fwrite(fd,zeros(1,30),'int32');
% c(1:3)

% the next two are supposed to be character strings; here the time is
% saved.
dateof = date;
dateout = [dateof(1:6),dateof(10:11)];
cnt=fwrite(fd,dateout,'char');
% d
% char(d(1:3))
% char(d(5:8))

% Two more ints...
cnt=fwrite(fd,[1,10],'int32');

% 10 more lines ....
%ns=min(e(2),10);
for i=1:10
	cnt=fwrite(fd,['#',char(zeros(1,78)),'#'],'char');
end;

% Now write the data.
cnt =fwrite(fd,vol,'float32');
fclose(fd);
fwt=1;
