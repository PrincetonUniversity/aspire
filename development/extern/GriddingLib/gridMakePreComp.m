function precomp=gridMakePreComp(n, nw)
% function precomp=gridMakePreComp(n, kernelwidth);
% Compute the np x 1 pre-compensation function, given the kernel width nw.
% If we are operating on an N-D array m, then np is obtained as:
%   sz=size(m); np=padfactor*sz(1);
%
% The basic idea of gridding is that the original (let's call it
% real-space) function is pre-compensated so that the "passband" errors of
% the interpolation kernel are cancelled.  To perform the compensation, the
% function computed here is multiplied by the real-space function along each dimension.
% The function is symmetric about np/2+1.

np=n*gridPadFactor('grid');
alpha=gridGetAlpha(nw);

if alpha>0

    % First, make an oversampled, frequency-space Kaiser function.
    ovc=4;  % oversamping
    kc=(-nw/2:1/ovc:nw/2-1/ovc)';  % nw*ovc points
    wc=kaiser(nw/2,alpha,kc);

    % Put the kaiser function into the middle of an np*ovc point array.
    spc=(np-nw)*ovc/2;  % Shift of origin for padding the array to np*ovc points.
    wcp=zeros(np*ovc,1);
    wcp(1+spc:nw*ovc+spc)=wc; % copy the Kaiser kernel.

    % Compute the FT (=IFT) and normalize it.
    fc=real(fftshift(fft(fftshift(wcp))));  % shift the center to the origin, ft and shift it back to the center.

    % Normalize the IFT and prepare to extract the central 1/ovc of the array.
    fc0=fc(np*ovc/2+1);   % pick up the zero-frequency point for normalization.
    spc=np*(ovc-1)/2;  % shift of origin to bring it to an array of size np

    % The compensation function is the reciprocal of the extracted IFT.
    precomp=fc0./fc(1+spc:np+spc);

else
    precomp=1;
end;
