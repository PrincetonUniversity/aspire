function mat=mrc2mat(mrcname)

stackreader=imagestackReader(mrcname);
n=stackreader.dim(3);
mat=stackreader.getImage(1:n);
