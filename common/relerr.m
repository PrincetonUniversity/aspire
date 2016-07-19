function err=relerr(A,B)
err= sqrt(sum((A(:)-B(:)).^2)/sum(A(:).^2));