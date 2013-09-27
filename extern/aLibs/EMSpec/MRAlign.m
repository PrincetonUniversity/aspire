function xytfr=MRAlign(imgs,refs,mode,tstep,doFlip,prior,sym,debug)
% function [xytfr alignedimages]=MRAlign(imgs,refs,mode,tstep, prior,sym);
% Multi-reference alignment of an image stack imgs to a stack of references
% refs.  It is up to the user to pre-mask the refs and images, so that they
% are zero at the periphery.
% The returned xytfr matrix is nimages x 6 and gives the xshift, yshift,
% theta, flip (0 or 1) and reference number.  The 6th element is the maximum
% cc value.  theta is in radians.  Mode =1
% for translate only, 2 for rotate only, 3 (default) for both.  tstep
% (default is 2) is the number of pixels at 80% radius for the size of the
% theta step (typically 2-4, but can be larger for speed).  Only part of
% the range of theta angles is searched if sym>1 (default is 1).
% At present the searches are not interpolated, and the shifts are
% quantized to 1 pixel.
% defaults are: mode=3, tstep=2, doFlip=1, sym=1, debug=0.

intShifts=0;

[n, ny, nim]=size(imgs);
[nxr, nyr, nrefs]=size(refs);

if nxr~=n  % change the refs size to match the images
    refs=Downsample(refs,n,1);
end;

if nargin<8
    debug=0;
end;
if nargin<7
    sym=1;
end;
if nargin<6
    prior=0;
elseif numel(prior)<2 && prior>0 % just a S.D.
    %     prior should be s_n^2 * |r-r0|^2/(2*s_r^2), where s_n is noise sigma,
    %     and s_r is translation sigma.  We assume s_n=1.
    prior=-Radius(n).^2/(2*prior^2);
end;
if nargin < 5
    doFlip=1;
end;
nflip=doFlip+1;

if nargin<4
    tstep=2;
end;
if nargin<3
    mode=3;
end;


if mode==1  % translate only
    ntheta=1;
    thetas=0;
    dt=0;
else
    if numel(tstep)==1
        ntheta=round(0.8*n*pi/(sym*tstep));  % rotation quantum is 2 pixels at 0.8 x R
        thmax=2*pi/sym;
        dt=thmax/ntheta;
        thetas=dt*((1-ntheta)/2:(ntheta-1)/2)';
    else
        ntheta=numel(tstep);
        thetas=tstep(:);
    end;
end;

if debug
    nrefs
    ntheta
    disp('Rotating references...');
end;
reflib=single(zeros(n,n,ntheta,nflip,nrefs));

for i=1:nrefs
    ref=refs(:,:,i);
    ref=ref/sqrt(ref(:)'*ref(:));  % normalize each reference
    for j=1:ntheta
        theta=thetas(j);
        if theta==0
            rref=ref;
        else
            rref=grotate(ref,theta);
        end;
        reflib(:,:,j,1,i)=rref;
        if nflip>1
            % flipud flips the first (x) coordinate, but doesn't move the fftcenter to
            % the right place if nx is even.  We flip the rotated reference.
            flipref=circshift(flipud(rref),[mod(n+1,2) 0]);
            reflib(:,:,j,2,i)=flipref;
        end;
    end;
end;

if mode ~=2
    freflib=single(complex(zeros(n,n,ntheta,nflip,nrefs)));
    for i=1:nrefs
        for j=1:ntheta;
            for l=1:nflip
                freflib(:,:,j,l,i)=conj(fftn(reflib(:,:,j,l,i)));
            end;
        end;
    end;
end;


if debug
    disp('Aligning...');
end;

% Do the alignment
% xytfr=zeros(nim,5);

xytfr=zeros(nim,6);

if mode==2  % rotate only
    cc=zeros(ntheta,1);
    mx=zeros(nflip,1);
    jx=zeros(nflip,1);
    mxl=zeros(nrefs,1);
    jt=zeros(nrefs,1);
    ifl=zeros(nrefs,1);
    for i=1:nim
        im=imgs(:,:,i);
        for j=1:nrefs
            for l=1:nflip
                for k=1:ntheta
                    rf=reflib(:,:,k,l,j);
                    cc(k)=im(:)'*rf(:);
                end;
                [mx(l), jx(l)]=max1di(cc);  % best theta for each ref
            end;
            [mxl(j), ifl(j)]=max(mx);  % Get the best over mirror images
            jt(j)=jx(ifl(j));         % Get the best rotation
        end;
        [mxcc, iref]=max(mxl);  % best ref overall
        xytfr(i,:)=[0 0 jt(iref) ifl(iref) iref mxcc];
        if mod(i,1000)==0
            disp(i);
        end;
    end;  % loop over images i
    
else
    
    % shift or shift-and-rotate.  The former simply has ntheta=1.
    ct=n/2+1;
    if debug
        disp('Rot and shift determination...');
    end;
    for i=1:nim
        cc1=single(zeros(n,n,ntheta));  % translational cc
        iflip=zeros(nrefs,1);
        coords=zeros(3,nrefs);
        mx=zeros(nflip,1);
        cmx=zeros(3,nflip);
        mxvals=zeros(nrefs,1);
        im=imgs(:,:,i);
        fim=fftn(im);
        for j=1:nrefs
            for l=1:nflip
                for k=1:ntheta
                    cc1(:,:,k)=fftshift(real(ifftn(freflib(:,:,k,l,j).*fim)))+prior;
                end;
                if ntheta>1
                    [mx(l), cmx(:,l)]=max3di(cc1);
                else
                    [mx(l), cmx(1,l), cmx(2,l)]=max2di(cc1);
                end;
            end;
            [mxvals(j), iflip(j)]=max(mx);
            coords(:,j)=cmx(:,iflip(j));  % encoded x,y,t for ref. j
        end; % loop over references j
        [mxcc, refj]=max(mxvals);
        theta=dt*(coords(3,refj)-(ntheta+1)/2);
        % q=[coords(1:2,refj)'-ct theta iflip-1 refj mxcc]
        xytfr(i,:)=[ct-coords(1:2,refj)' theta iflip(refj)-1 refj mxcc];
        
        if mod(i,1000)==0
            disp(i);
        end;
    end;  % loop over images i
end;

%
% % if nargout>1  % we return aligned images
% if debug
%     disp('Rotating output images...');
% end;
% alignedimages=zeros(n,n,nim);
% for j=1:nim
%     im=imgs(:,:,j);
%     theta=-xytfr(j,3);
%     iflip=xytfr(j,4);
%     if intShifts
%         im=circshift(imgs(:,:,j),-round(xytfr(j,1:2)));
%     else
%         im=shiftf(imgs(:,:,j),-xytfr(j,1:2));
%     end;
%     if theta ~= 0
%         im=grotate(im,theta);
%     end;
%     if iflip
%         im=circshift(flipud(im),[mod(n+1,2) 0]);  % flip after rotate and shift.
%     end;
%     alignedimages(:,:,j)=im;
% end;
%
