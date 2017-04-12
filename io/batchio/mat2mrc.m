function mat2mrc(mrcname,mat)

n=size(mat,3);
stackwriter=imagestackWriter(mrcname,n);
stackwriter.append(mat);
stackwriter.close;