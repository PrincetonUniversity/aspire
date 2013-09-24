function p1=gridExtractLine(P, theta, kernelsize);
% function p1=gridExtractLine(P, theta, kernelsize);
% Compute the central-section line of the ft of an image,
% using the data structure P from the function gridMakePaddedFT.
% The line is returned as the data structure p1.  The extracted line is in
% p1.PadFT and has p1.np=P.np points and is returned
% in cosine-and-sine form, i.e. a real vector.
% Theta is the angle ccw from the x-axis.
% Line extraction takes 3 ms for n=64 points, 5 ms for 128 points
% and kernelsize=5 on my 1.5GHz powerbook.
% Vectorized version.  fs 18 Oct 06

% parameters and static variables
ov=1024;  % oversampling of window for look-up
persistent w1 % Cache for the w1 lookup table.

% default argument
if nargin<3
    kernelsize=5;
end;

nw=kernelsize;
alpha=gridGetAlpha(nw);
if alpha==0
    error('Invalid kernel size; it must be odd.');
end;


if (numel(w1)~=nw*ov) % if the cached table is absent or not the right size
    % Create the interpolation kernel
    dw=1/ov;  % fractional increment.  We assume ov is even
    k=(-nw/2+dw/2:dw:nw/2-dw/2)';  % nw*ov space points
    w=kaiser(nw/2,alpha,k);  % Compute the Kaiser window
    w=w*ov/sum(w);  % normalize it.
    w1=zeros(nw,ov);

    % Make the 1D lookup table w1(:,i).  i=1 corresponds to a shift between -0.5 and
    % 0.5+dw, so we assign it the mean shift -0.5+dw/2.
    % i=ov corresponds to a mean shift of 0.5-dw/2.
    for i=1:ov
        w1(:,i)=w(ov-i+1:ov:nw*ov);
    end;
end;


% Interpolation is done here.

np=P.np;
np1=P.np1;
np1ctr=P.np1ctr;

% % Old non-vectorized code is here:
% rot=[c s;-s c];  % rotation matrix
% sp1=P.sp1;
% iline=zeros(np1,1);
% nw2=floor(nw/2);  % half-width of kernel (=(nw-1)/2).
% i=np1ctr;  % Compute the central section, rotated from the y-axis
% for j=np1ctr-np/2:np1ctr+np/2-1  % get np points
%     p=rot*([i j]'- centerp)+centerp;
%     pint=round(p);
%     pfrac=floor((p-pint)*ov)+ov/2+1;
%     i1=pint(1);
%     j1=pint(2);
%     iline(j)=w1(:,pfrac(1))'*P.PadFT(i1-nw2:i1+nw2,j1-nw2:j1+nw2)*w1(:,pfrac(2));
% end;
% line=iline(sp1+1:sp1+np);

% Vectorized code.
p1=gridMakeNullFT(P.n,1);
pft=zeros(np,1);
nw2=floor(nw/2);  % half-width of kernel (=(nw-1)/2).
% Compute the central section.  The j's are implicitly set to 0, to rotate
% from the x-axis.
js=(-np/2:np/2-1)';  % we compute np points along the y' axis.  Always, x (called is)=0.
ip=sin(theta)*js+np1ctr;
jp=cos(theta)*js+np1ctr;
ipint=round(ip);
jpint=round(jp);
ipfrac=floor((ip-ipint)*ov)+ov/2+1;
jpfrac=floor((jp-jpint)*ov)+ov/2+1;
addrs=ipint+np1*(jpint-1);  % 1-dim addresses in the PadFT array.
for k=-nw2:nw2
    for l=-nw2:nw2
        pft=pft+w1(k+nw2+1,ipfrac)'.*w1(l+nw2+1,jpfrac)'.*P.PadFT(addrs+k+np1*l);
    end;
end;
sp1=p1.sp1;
% Return the padded ft with np1 points.
p1.PadFT=zeros(np1,1);
p1.PadFT(1+sp1:np+sp1)=pft;
