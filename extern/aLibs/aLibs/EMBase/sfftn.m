function F=sfftn(f)
F=fftshift(fftn(ifftshift(f)));
