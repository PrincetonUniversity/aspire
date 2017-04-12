function mat=mrc2mat_complex(mrcname)

stackreader=imagestackReaderComplex(mrcname);
n=stackreader.dim(3);
mat=stackreader.getImage(1:n);