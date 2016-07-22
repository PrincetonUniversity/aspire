function MakeJpegsFromEMFiles(OutputDir, binning, display)
% function MakeJpegsFromEMFiles(binning, display)
% For each EM image file (.mrc, .dm3, .img, .hed, .tif)
% in the current directory, make a .jpg file.  All arguments are optional.
% If OutputDir is given and is not an empty string, the
% jpeg files are written there.  The argument binning controls the degree
% of downsampling.  Its default value is 2.  If display=0 then no display
% is shown.
% fs July 2009 rev. Apr 2011

disp(['Converting EM files to jpgs in directory ' pwd]);

if nargin<1
    OutputDir='';
end;
len=numel(OutputDir);
if len>0
    if OutputDir(len)~='/' % doesn't end with slash
        OutputDir=[OutputDir '/'];
    end;
    disp(['Writing output files to ' OutputDir]);
end;
if nargin<2
    binning=2;
end;
if nargin<3
    display=1;
end;

d=dir;
for i=3:numel(d)
    [mi pixA ok]=ReadEMFile(d(i).name);
    if ok
        disp(d(i).name);
        m=RemoveOutliers(mi,4);
        n=size(m,1);
        if binning>1
            m=Downsample(m,n/binning);
        end;
        ms=uint8(imscale(m));
        if display
            figure(1); SetGrayscale;
            imac(ms);
            title([d(i).name '  pixA= ' num2str(pixA*binning)]);
            drawnow;
        end;
        [pa nm]=fileparts(d(i).name);
        imwrite(ms,[OutputDir nm '.jpg'],'jpg');
    else
        skip=d(i).name
    end;
end;
disp('Done.');
