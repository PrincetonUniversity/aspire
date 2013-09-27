function MakeJpegsFromDm3s(binning)
% function MakeJpegsFromDm3s(binning)
% For each .dm3 file in the current directory, make a .jpg file.  The
% optional argument binning controls the degree of downsampling.  Its
% default value is 2.
% fs July 2009

disp(['Converting dm3s to jpgs in directory ' pwd]);
if nargin<1
    binning=2
end;
d=dir;
for i=3:numel(d)
    p=strfind(d(i).name,'.dm3');
    if numel(p)>0
        basename=d(i).name(1:p-1);
        disp(d(i).name);
        m=single(ReadDM3(d(i).name));
        m=RemoveOutliers(m,4);
        n=size(m,1);
        if binning>1
            m=Downsample(m,n/binning);
        end;
        ms=uint8(imscale(m));
        imwrite(ms,[basename '.jpg'],'jpg');
    end;
end;
disp('Done.');
