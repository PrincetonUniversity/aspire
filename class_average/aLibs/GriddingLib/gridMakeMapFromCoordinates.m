% function rvol=gridMakeMapFromCoordinates(n, coords, weights, sigma)
% function rvol=gridMapFromCoordinates(n, coords, weights, sigma)
% Very fast creation of 3D map from pdb-style coordinates.  The map is
% created in an n x n x n volume; each atom is represented by a Gaussian
% density of standard deviation sigma.  The coordinates are a 3 x na
% matrix of coordinates in A. The weights (na x 1) are the amplitudes
% given to each atom.
% If sigma is zero, a high-resolution map (3A) is created using a
% sharp-cutoff filter.
% n must be a multiple of 8, and to avoid truncation artifacts, it's best
% that coordinate values are between 2 and n-1.

% Test code:
n=128;
na=2;
% coords=2+(n-3)*rand(3,na);
for jj=0:.1:1
    coords=ones(3,2)*n/2+1;
    coords(1,1)=n/2+1+jj;
    coords(1,2)=n/2-8;
    weights=ones(na,1);
    sigma=-1;
    % disp([num2str(na) ' atoms']);
    % tic
    % R3=gridMakeMapFromCoordinates(n, coords, weights, sigma);
    % toc
    % figure(1); SetGrayscale;
    % imacs(sum(R3,3));  % display the volume summed along z
    %
% Actual function
    % parameters
    na=size(coords,2);
    mode='grid';
    kernelsize=5;
    gridfactor=1.25;
    
    nw=kernelsize;  % size of interp window
    nw2=floor(nw/2);  % half-width of kernel, =1 for nw=3.
    
    np=n*gridfactor;  % padded size
    np1=np+NextNiceNumber(2*nw2,7,4);  % expanded to include a border for kernel points
    
    % Make the look-up table for interpolation
    [w1 ov]=gridMakeKaiserTable(kernelsize,mode);
    
    % We use double-precision everywhere to avoid artifacts in the FT
    % computation of Gaussians.
    pvol=zeros(np1,np1,np1);  % Real-space padded volume
    
    sp1=(np1-np)/2;  % offset into np1-sized volume.
    
    ovctr=ov/2+1;  % Center of oversampled interpolation table.
    cube=zeros(nw,nw,nw);
    
    % We expand the coordinates to fill an np-sized volume, with sp1 points of
    % padding.
    coords=(coords-1)*gridfactor+sp1+1;
    cint=round(coords);
    cfrac=floor((coords-cint)*ov)+ovctr;
    
    for ia=1:na
        % volume coordinates
        i1=cint(1,ia);
        j1=cint(2,ia);
        k1=cint(3,ia);
        % we compute the interpolated function on a cube surrounding (i1 j1 k1)
        % one point at a time.
        insplane=weights(ia)*w1(:,cfrac(1,ia))*w1(:,cfrac(2,ia))';
        for k2=1:nw
            cube(:,:,k2)=insplane*w1(k2,cfrac(3,ia));
        end;
        % Add the cube into the volume.
        pvol(i1-nw2:i1+nw2,j1-nw2:j1+nw2,k1-nw2:k1+nw2)...
            =pvol(i1-nw2:i1+nw2,j1-nw2:j1+nw2,k1-nw2:k1+nw2)+cube;
    end; %for ia
    
    pvol=Crop(pvol,np);  % get rid of the border.
    fvol=fftshift(fftn(pvol));  % zero-centered FT
    
    % We will filter the volume by both the post-compensation function
    % (high-frequency boost) and a gaussian or (high resolution) fuzzymask.
    postcomp=gridMakePreComp(n,nw);  % makes an np x 1 vector, zero centered
    pc2=kron(postcomp,postcomp');
    pc3=reshape(kron(pc2(:),postcomp'),[np np np]);  % 3D compensation function
    if sigma>0
        pc3=pc3.*exp(-2*(pi*sigma*double(Radius3(np))/np).^2 ); % multiply by Gaussian
    elseif sigma<0  % exponential function
        pc3=pc3./(1+(2*pi*sigma*double(Radius3(np)/np)).^2);
        %     pc3=pc3.*(fuzzymask(np,3,.35*n,.2*n));  % half power at sampling freqency/3
        sigma0=1.2;
        pc3=pc3.*exp(-2*(pi*sigma0*double(Radius3(np))/np).^2 ); % multiply by Gaussian
        pc3=pc3/sum(pc3(:))*n^3;
    else  % sigma=0: sharp filter.
        pc3=pc3.*sqrt(fuzzymask(np,3,n/3,.05*n));  % half power at sampling freqency/3
    end;
    
    rvol=real(ifftn(fftshift(Crop(fvol.*pc3,n))));
    
    subplot(2,2,1);
    imacs(((sum(rvol,3))));
    subplot(2,2,3);
    % semilogy(-n/8:n/8-1,Crop(sect(rvol),n/4),'k.-');
    plot(-n/8:n/8-1,log(Crop(sect(rvol),n/4)),'k.-');
    drawnow;
end;
return
%
% q=sect(rvol);
% semilogy(1:n/2,abs(Crop(q,n/2)),'k.-');
% axis([-inf inf 1e-9 1e-1]);
% % plot(1:n/2,Crop(q,n/2),'k.-');
% % axis([-inf inf -3e-3 3e-3]);

