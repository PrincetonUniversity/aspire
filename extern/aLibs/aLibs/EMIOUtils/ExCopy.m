function ExCopy(infile, outfile, ptrs)
% Exclusive copy of Imagic images.
% Let j=ptrs(i). For each i, copy the jth image from infile to the ith
% position in outfile.
% We make assumptions: only 1 header record per image, and all images are
% the same size.
% The ptrs vector can be gotten thus: ptrs=find(1-exclude), where exclude
% is a vector that is 0 for ok, 1 for excluding a given image in the input
% file.

[inhed inimg]=makeImagicFileNames(infile);
[outhed outimg]=makeImagicFileNames(outfile);

% Read the input header

% disp('Reading the input header')

IHdr=ReadImagicHeader(inhed);

iimg=fopen(inimg,'r',IHdr.ByteOrder);  % Open the image data
oimg=fopen(outimg,'w',IHdr.ByteOrder);

ni=numel(ptrs);
% ni

OHdr.Vals=IHdr.Vals(:,1:ni);  % Pre-allocate the output.
OHdr.Name(1:ni,:)=IHdr.Name(1:ni,:);
OHdr.Strings=IHdr.Strings(1:ni,:);
OHdr.Type=IHdr.Type(1,:);
OHdr.ByteOrder=IHdr.ByteOrder;

% disp('Copying headers and images');

    OHdr.Type=IHdr.Type;
    [bytes elems type]=PixelSize(OHdr.Type);

    for i=1:ni
    j=ptrs(i);  % pointer to the input image.
    OHdr.Vals(:,i)=IHdr.Vals(:,j);  % Copy the header information over.
    OHdr.Name(i,:)=IHdr.Name(j,:);
    OHdr.Strings(i,:)=IHdr.Strings(j,:);
    
    % correct some header entries
    OHdr.Vals(1,i)=i;
    OHdr.Vals(2,i)=ni-i;

    % figure out the size of the image data


%     if i==1; type,bytes,elems, end

    imageBytes=bytes*OHdr.Vals(12,i);
    imageElems=elems*OHdr.Vals(12,i);

    ok=fseek(iimg,imageBytes*(j-1),'bof'); % seek to image data
    imdata=fread(iimg,imageElems,['*' type]);
    count=fwrite(oimg,imdata,type);
%     if i<10; count, end
end;

fclose(iimg);
fclose(oimg);
WriteImagicHeader(OHdr,outhed);

end



function [hdrname imgname]=makeImagicFileNames(basename)
% Ignore the extension, and construct the filenames *.hed and *.img
[n1 slen]=size(basename);
if (slen>4) && (basename(slen-3)=='.') % has an extension
    %     if (basename(slen-2:slen)=='hed') || (basename(slen-2:slen)=='img')
    if strcmp(basename(slen-2:slen),'hed') || strcmp(basename(slen-2:slen),'img')
        basename=basename(1:slen-4);  % remove the extension.
    end;
end;
hdrname=strcat(basename,'.hed');
imgname=strcat(basename,'.img');
end


function [bytes elems type]=PixelSize(typeString)
elems=1;
switch typeString
    case 'PACK'
        bytes=1;
        type='uint8';
    case 'INTG'
        bytes=2;
        type='int16';
    case 'REAL'
        bytes=4;
        type='float32';
    case 'COMP'
        bytes=8;
        type='float32';
        elems=2;
    case 'RECO'
        bytes=8;
        type='float32';
        elems=2;
    otherwise
        error(['ExCopy: unsupported Imagic data mode: ',typeString]);
        bytes=0;
        type='???';
        elems=0;
end
end
