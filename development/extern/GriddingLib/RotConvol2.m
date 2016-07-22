function mr=RotConvol2(m,rotmat,table,radius,inctr,rctr)
% function mr=RotConvol2(m,rotmat,table,radius,inctr,rctr)
% Convolution loop used
% in grotate.  Given the input image m, use the inverse rotation matrix rotmat and
% the kernel look-up table (nw x ntcols in size) to perform the convolution
% by the nw x nw kernel (nw must be odd).  Only a disc region of size
% radius in the output image mr will be set (the rest are zeros).  The
% center locations in the input image and the corresponding point in the
% output image are inctr and rctr, respectively.  No array bound checks are made, so
% the size of m must be sufficient so that inctr and rctr must be a
% distance of at least radius + nw2 from each edge of the image, where nw2
% is the kernel half-width.

[ni ni2]=size(m);
mr=zeros(ni,ni2);
% [nr nr2]=size(mr);
[nw ntcols]=size(table);
nw2=(nw-1)/2;   % kernel is assumed to have odd width.
% ntrows is equal to the table oversampling factor.

r0=trunc(radius);
for dj=-r0:r0  % dj is the distance in the second dimension from the center (input)
    r1=trunc(sqrt(radius^2-dj^2));
    ir=rctr(1)+rotmat(1,2)*dj-rotmat(1,1)*r1;  % 1st dimension of output
    jr=rctr(2)+rotmat(2,2)*dj-rotmat(1,2)*r1;  % 2nd dim of output
    j=inctr(2)+dj;  % absolute 2nd dim (input)
    dri=rotmat(1,1); % shift of output per 1st dim step of input.
    drj=rotmat(2,1);
    
    for di=-r1:r1  % loop over 1st dimension of input
        i=inctr(1)+di;
        iri=round(ir); % pick up integer and fraction parts.
        irf=nw*ntcols*floor(ir-iri+0.5)+nw2+1;  % offset so we can add -nw2 to get first element.
        jri=round(jr);
        jrf=nw*ntcols*floor(jr-jri+0.5)+nw2+1;
% note that ir and jr are both positive, so rounding isn't too difficult.
        sum2=0;
        for jk=-nw2:nws % Convolution loop
            sum1=0;
            for ik=-nw2:nw2;
                sum1=sum1+m(iri,jri)*table(irf+ik);
            end;
            sum2=sum2+sum1*table(jrf+jk);
        end;
        mr(di+inctr(1),dj+inctr(2))=sum2;
        
        ir=ir+dri;
        jr=jr+drj;
    end;
end;
