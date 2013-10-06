function c=CrossCorrelation(a,b)
% Compute the cross-correlation function.
% For alignment: if the maximum of c is found to be
% at (i,j) then b=shift(b,i,j) will bring b into
% alignment with a.

a=a-mean(mean(a));
b=b-mean(mean(b));

% ja=size(a,2);
% jb=size(b,2);

% % Pad the arrays to make them the same size
% if (ia > ib)
% 	b(ia,jb)=0;	% Force b to have as many rows as a
% elseif (ib > ia)
% 	a(ib,ja)=0;
% end;
% if (ja > jb)
% 	b(ib,ja)=0;	% Force b to have as many columns as a
% elseif (jb>ja)
% 	a(ia,jb)=0;
% end;

fa=fft2(a);
fb=fft2(b);
c=real(ifft2(fa.*conj(fb)));
