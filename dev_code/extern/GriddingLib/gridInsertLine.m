function P1=gridInsertLine(p1, P, theta, kernelsize);
% function P1=gridInsertLine(p1, P, theta, kernelsize);
% Add the central-section line to the ft of an image,
% using the data structure P from the function gridMakePaddedFT.
% The line is oriented at the angle theta ccw from the x axis.
% The line is assumed to be
% in cosine-and-sine form, i.e. a real vector.
% Vectorized code.

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

P1=P;
np=P.np;
np1=P.np1;
np1ctr=P.np1ctr;
sp1=P.sp1;

% Old scalar code:
s=sin(theta);
c=cos(theta);
rot=[c s;-s c];  % rotation matrix
iline=p1.PadFT;
% 
centerp=[np1ctr np1ctr]';  % coordinate of the origin
nw2=floor(nw/2);  % half-width of kernel (=(nw-1)/2).
i=np1ctr;  % Compute the central section
for j=np1ctr-np/2:np1ctr+np/2-1  % get np points
    p=rot*([i j]'- centerp)+centerp;
    pint=round(p);
    pfrac=floor((p-pint)*ov)+ov/2+1;
    i1=pint(1);
    j1=pint(2);
    P1.PadFT(i1-nw2:i1+nw2,j1-nw2:j1+nw2)=P1.PadFT(i1-nw2:i1+nw2,j1-nw2:j1+nw2)...
        +iline(j)*w1(:,pfrac(1))*w1(:,pfrac(2))';
end;



% % Vectorized code.  It doesn't work.
% % 
% nw2=floor(nw/2);  % half-width of kernel (=(nw-1)/2).
% % Compute the central section.  The j's are implicitly set to 0.
% is=(-np/2:np/2-1)';  % we will insert np points, one for each element of line.
% ip=cos(theta)*is+np1ctr;  % x-coordinate of each point in the PadFT
% jp=sin(theta)*is+np1ctr;  % y-coordinate of each point.
% ipint=round(ip);
% jpint=round(jp);
% ipfrac=floor((ip-ipint)*ov)+ov/2+1;  % Fractional part is multiplied by ov for lookup.
% jpfrac=floor((jp-jpint)*ov)+ov/2+1;
% addrs=ipint+np1*(jpint-1);  % Construct the 1-dim addresses for the PadFT array.
% for k=-nw2:nw2  % loop over each kernel element
%     for l=-nw2:nw2
%         addrp=addrs+k+np1*l;  % modify the addresses to do the integer shifts by k and l.
%         
%          Accum(line.*w1(k+nw2+1,ipfrac)'.*w1(l+nw2+1,jpfrac)',P.PadFT,addrp);
% %         P1.PadFT(addrp)=P.PadFT(addrp)+line.*w1(k+nw2+1,ipfrac)'.*w1(l+nw2+1,jpfrac)';
% % for comparison, here is the extract-line code.
%         %         line=line+w1(k+nw2+1,ipfrac)'.*w1(l+nw2+1,jpfrac)'.*P.PadFT(addrs+k+np1*l);
%     end;
% end;
