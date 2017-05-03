function mat2mrc_complex(mrcname,mat)

n=size(mat,3);
stackwriter=imagestackWriterComplex(mrcname,n);
stackwriter.append(mat);
stackwriter.close;