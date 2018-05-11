function cryo_raynormalize_outofcore(inname,outname)
%
% CRYO_RAYNORMALIZE_OUTOFCORE  Normalize a Fourier rays to energy 1.
% 
% cryo_raynormalize_outofcore(inname,outname)
%   Normalize a dataset of Fourier rays so that each ray has energy 1.
%   inname is the name of a comllex-valued stack (generated, for example,
%   by cryo_pft_outofcore). outname is the name of the stack to store the
%   normalized data.
%
%   See cryo_raynormalized.m for more information.
%
% Yoel Shkolnisky, April 2017.
 
in=imagestackReaderComplex(inname);
n_projs=in.dim(3);
n_theta=in.dim(2);
out=imagestackWriterComplex(outname,n_projs);

for k=1:n_projs
    pf=in.getImage(k);
    for j=1:n_theta
        nr=norm(pf(:,j));
        if nr<1.0e-13
            warning('Ray norm is close to zero. k=%d  j=%d',k,j);
        end
        pf(:,j)=pf(:,j)./nr;
    end
    out.append(pf);
end
out.close;