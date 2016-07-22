function cc=CTFCrossCorr(spect,res,limits,P,sqroot)
% Compute the cross correlation of the centered spectrum spect with the CTF
% function computed from the CTF parameters P.

if nargin<5
    sqroot=0;
end;

nu=size(spect,1);
if sqroot
    ct=abs(CTF(nu,res,P)).*limits;
else
    ct=CTF(nu,res,P).^2.*limits;
end;
ct=ct(:);
spectm=limits.*spect;
spectm=spectm(:);

cc=ct'*spectm/sqrt((ct'*ct)*(spectm'*spectm));
