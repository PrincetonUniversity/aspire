function fwt = writeMRC(vol, pixelsize, filename);
% function fwt = writeMRC(vol, pixelsize, filename): write a 3D map into the mrc
% format (little-endian assumed).
fwt=0;
% Here is where we specify little-endian data:
fd=fopen(filename,'wb','ieee-le');
[nx, ny, nz]=size(vol);

% the first 10 values, which are integers:
% nc nr ns mode ncstart nrstart nsstart nx ny nz
a(1:10)=[nx,ny,nz,1,nx/2+1,ny/2+1,nz/2+1, nx,ny, nz];
cnt=fwrite(fd,a,'int32');
% a(1:10)

% the next 12, which are floats.
% the first three are the cell dimensions.
% xlength ylength zlength alpha beta gamma mapc mapr maps
% amin amax amean.
b=[nx*pixelsize,ny*pixelsize,nz*pixelsize,zeros(1,9)];
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
