% SplitMRCFiles.m
%  Select mrc files from SerialEM and split them into individual mrc
%  images.  Also writes out jpeg images binned by 2.

ndigits=3; % number of digits appended
separator='_';  % could be '_' for example.
WriteJpegs=1;

figure(1);
SetGrayscale;
uiPath=pwd;
[files uiPath]=uigetfile('*','Select mrc files','multiselect','on',uiPath);
if isnumeric(files);  % returns 0 if cancelled
    return
end;
if ischar(files)
    inputFilenames{1}=files;
    nim=1;
else
    inputFilenames=files;
    nim=numel(files);
end;
cd(uiPath);

% uiPath
disp('Input files');
disp(char(inputFilenames));
micrographPath='Micrograph/';
jpegPath='Micrograph/Jpeg/';
if ~DirectoryExists([uiPath micrographPath])
    mkdir([uiPath micrographPath]);
end;
if ~DirectoryExists([uiPath jpegPath])
    disp(['Making the directory ' uiPath jpegPath]);
    mkdir([uiPath jpegPath]);
end;


for i=1:nim
    name=inputFilenames{i};
    fullName=[uiPath name];
    [m s]=ReadMRC(fullName);
    if isa(m,'uint8')
        mode=0;
    elseif isa(m,'int16')
        mode=1;
    elseif isa(m,'single')
        mode=2;
    elseif isa(m,'uint16')
        mode=6;
    else
        mode=2;  % default is float
    end;
    nparts=size(m,3);
    [pa baseName ex]=fileparts(name);
    if nparts<2
        disp('This file contains a single image');
    end;
    disp('Writing files:')
    for j=1:nparts
        endString=sprintf([separator '%0*d'],ndigits,j);
        outName=[uiPath micrographPath baseName endString '.mrc'];
        disp(outName);
        WriteMRC(m(:,:,j),s.pixA,outName,mode);
        %                     outName=[uiPath baseName endString '.mrc'];
        
        mj=BinImage(m(:,:,j),2);
        imacs(mj);
        title([baseName endString],'interpreter','none');
        drawnow;
        if WriteJpegs
            jpegName=[uiPath jpegPath baseName endString '.jpg'];
            %            imwrite(uint8(imscale(mj)),jpegName);
            imwrite(uint8(imscale(mj,256,1e-3)),jpegName);
        end;
    end;
end;
