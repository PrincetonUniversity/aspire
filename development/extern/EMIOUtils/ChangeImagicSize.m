function byteorder=ChangeImagicSize(infile,outfile,n,HeaderOnly)
% function ChangeImagicSize(infile,outfile,n [,type [,HeaderOnly]])
% The output file is of type REAL.
% Change each (square) image in an Imagic image stack to be n x n in size,
% using the Downsample function to shrink or expand it.
% If HeaderOnly=1, only the header is changed; it is up to the user to
% modify the .img file to correspond.

    type='REAL';
if nargin<4
    HeaderOnly=0;
end;

% 
% n=64;
% infile='Repicked'
% outfile='shrink'

chunk=1024;

% images to operate on at a time.


[path base ext v]=fileparts(infile);
[patx baso exo v]=fileparts(outfile);
% base
% baso
inhedname=fullfile(path, [base '.hed']);
outhedname=fullfile(path, [baso '.hed']);
disp('Change Imagic Size');
disp(['Reading header ' inhedname]);
Hdr=ReadImagicHeader(inhedname);

byteorder=Hdr.ByteOrder;
OrigType=Hdr.Type(1,:);
sizes=Hdr.Vals(12:14,1);
n1=sizes(2);
n0=sizes(3);
[nn nim]=size(Hdr.Vals);

for i=1:nim
    Hdr.Vals(12:14,i)=[n*n n n]';
    Hdr.Type(i,:)=type;
end;

disp(['Writing header ' outhedname]);
WriteImagicHeader(Hdr,outhedname);

if HeaderOnly
    return
end;

% Read and write the image data
disp(['Now, read and write ' num2str(nim) ' images']);
inimgname=fullfile(path, [base '.img'])
outimgname=fullfile(path, [baso '.img'])

inimg=fopen(inimgname,'r',byteorder);
outimg=fopen(outimgname,'w',byteorder);

nel=chunk*n0*n1;
done=0;
total=0;
if strcmp(OrigType,'PACK')
    format='uint8=>float32';
else
    format='*float32';
end;
while ~done
    images=fread(inimg,nel,format);
    nim=numel(images)/(n0*n1);
    images=reshape(images,n1,n0,nim);
    done = (nim < chunk);
    if nim==0 
        break
    end;
    
    dimages=Downsample(images,n,1);
    count=fwrite(outimg,dimages,'float32');
    total=total+count/(n*n);
    ImagesWritten=total
end;

fclose(inimg);
fclose(outimg);
disp('Done.');
