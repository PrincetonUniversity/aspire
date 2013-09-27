function out=BinImage(in,nb)
% function out=BinImage(in,nb)
% Bin an image (or a stack of them) by averaging nb x nb pixels to
% produce each output pixel.  If the image size is not a multiple of nb, it
% is cropped first.  Rectangular images are handled correctly.

if nb<=1
    out=in;
    return;
end;
nb=round(nb);

dims=ndims(in);
[nx ny nim]=size(in);
n=[nx ny];
dims=2;
if any(n==1)  % a 1d function
    dims=1;
    n=prod(n);
    in=reshape(in,n,nim);
end;
if numel(nb)<2  % if nb is a scalar, make it a vector
    nb=nb(1)*ones(1,2);
end;

% Check to see if n is a multiple of nb, and if not, force it.
if any(mod(n,nb))
    n=floor(n./nb).*nb;
    in=Crop(in,n);
end;

% Do the binning
n1=n./nb;
out=zeros([n1 nim]);
switch dims
    case 1
        for i=1:nim
            % Bin along x
            in1=reshape(in(:,i),nb,n1);   % make nb adjacent pixels a column
            out(:,i)=mean(in1)';
        end;
    case 2
        for i=1:nim
            % Bin along x
            nx=n1(1)*n(2);
            in1=reshape(in(:,:,i),nb(1),nx);   % make nb adjacent pixels a column
            out1=mean(in1);
            out1=reshape(out1,n1(1),n(2))';  % convert back to rectangular image
            % Bin along y
            ny=prod(n1);
            out1=reshape(out1,nb(2),ny); % take the mean along y
            out(:,:,i)=reshape(mean(out1),n1(2),n1(1))';
        end;
    otherwise
        error('Number of dimensions out of range (1-2)');
end;