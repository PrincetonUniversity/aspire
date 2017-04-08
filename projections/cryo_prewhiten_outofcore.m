function cryo_prewhiten_outofcore(instackname, outstackname,noise_response )
% CRYO_PREWHITEN_OUTOFCORE Pre-whiten a stack of images.
%  
% cryo_prewhiten_outofcore(instackname,outstackname,noise_response)
%   Prewhiten projections in the MRC file named instackname, and write the
%   prewhitened images to the MRC file outstackname. 
%   The images are prewhitened using the filter noise_response.
%   The function does not load the entire stack into memory.
%
% Example:
%   cryo_prewhiten_outofcore('instack.mrc','outstack.mrc',noise_response);

% See cryo_prewhiten for more details.
%
% Yoel Shkolnisky, May 2016.

delta=eps('single');

instack=imagestackReader(instackname);
n=instack.dim(3);
L=instack.dim(1);
l=floor(L/2);
K=size(noise_response, 1);
k=ceil(K/2);

% The whitening filter is the sqrt of of the power spectrum of the noise.
% Also, normalized the enetgy of the filter to one.
filter=sqrt(noise_response);      
filter=filter/norm(filter(:));

% The power spectrum of the noise must be positive, and then, the values
% in filter are all real. If they are not, this means that noise_response
% had negative values so abort.
assert(norm(imag(filter(:)))<10*delta); % Allow loosing one digit.
filter=real(filter);  % Get rid of tiny imaginary components, if any.

% The filter should be cicularly symmetric. In particular, it is up-down
% and left-right symmetric.
assert(norm(filter-flipud(filter))<10*delta); 
assert(norm(filter-fliplr(filter))<10*delta);

% Get rid of any tiny asymmetries in the filter.
filter=(filter+flipud(filter))./2;
filter=(filter+fliplr(filter))./2;

% The filter may have very small values or even zeros. We don't want to
% process these so make a list of all large entries.
nzidx=find(filter>100*delta);
fnz=filter(nzidx);

outstack=imagestackWriter(outstackname,n,1);

%fprintf('Whitening...\n');
printProgressBarHeader;
for i=1:n
    progressTic(i,n);
    pp=zeros(K);
    proj=instack.getImage(i);
    if mod(L,2)==1 % Odd-sized image
        pp(k-l:k+l, k-l:k+l)=proj; % Zero pad the image to twice the size.
    else
        pp(k-l:k+l-1, k-l:k+l-1)=proj; % Zero pad the image to twice the size.
    end
    fp=cfft2(pp); % Take the Fourier transform of the padded image.
    p=zeros(size(fp));
    p(nzidx) = fp(nzidx)./fnz; % Divide the image by the whitening filter, 
                               % but onlyin places where the filter is
                               % large. In frequnecies where the filter is
                               % tiny  we cannot pre-whiten so we just put
                               % zero.
    p2 = icfft2(p);
    assert(norm(imag(p2(:)))/norm(p2(:))<1.0e-13); % The resulting image should be real.
    if mod(L,2)==1
        p2 = p2(k-l:k+l, k-l:k+l);   
    else
        p2 = p2(k-l:k+l-1, k-l:k+l-1);
    end
    outstack.append(real(p2));
end;

outstack.close
end

