function c=FromCAS(r)
% Convert the "cosine and sine" real array r back to complex c.  We assume
% the origin is at the center of the arrays (i.e. fftshift has been
% applied.)  The -n/2 point is set to zero, which is at c(1,1).
% This routine handles arrays up to 3D

[s1, s2, s3]=size(r);
o1=min(2,s1);
o2=min(2,s2);
o3=min(2,s3);
c=zeros(s1,s2,s3);
c(o1:s1,o2:s2,o3:s3)=0.5*(r(o1:s1,o2:s2,o3:s3)+r(s1:-1:o1,s2:-1:o2,s3:-1:o3)...
    +1i*(r(o1:s1,o2:s2,o3:s3)-r(s1:-1:o1,s2:-1:o2,s3:-1:o3)));
