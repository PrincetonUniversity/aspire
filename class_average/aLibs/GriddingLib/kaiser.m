function w=kaiser(n,alpha,k)
% function w=kaiser(n,alpha,k)
% Kaiser-Bessel window function as defined in Hamming's book.
% n is the half-width of the final window; k is the list of window values.
r0=besseli(0,alpha);
w=k*0;
for i=1:numel(k)
  q=k(i);
  if abs(q)>n
    w(i)=0;
  else
    r=sqrt(1-(q/n).^2);
    w(i)=besseli(0,alpha*r)./r0;
  end;
end;
