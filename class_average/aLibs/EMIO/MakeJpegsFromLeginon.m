% MakeJpegsFromLeginon
figure(1);
clf; SetGrayscale;
sourcepath=pwd;
disp(['Converting mrc images to jpgs from directory ' sourcepath]);
cd ../jpegs
destpath=pwd;
disp(['to ' destpath]);
cd(sourcepath);
decimation=2;
d=dir;
for i=3:numel(d)
    [path basename ext]=fileparts(d(i).name);
    if strcmpi(ext,'.mrc')
        disp(d(i).name);
        [m s]=ReadMRC(d(i).name);
        pixelsize=s.rez/s.nx;
        m=RemoveOutliers(single(m),4);
        n=size(m,1);
        if decimation>1
            m=BinImage(m,decimation);
        end;
        ms=uint8(imscale(m));
         image(ms);
        title(d(i).name,'interpreter','none');
        drawnow;
        imwrite(ms,[destpath '/' basename '.jpg'],'jpg');
    end;
end;
disp('Done.');
