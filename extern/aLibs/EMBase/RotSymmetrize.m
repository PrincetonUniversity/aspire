function out=RotSymmetrize(in,angBW)
% function out=RotSymmetrize(in, angBW);
% Rotationally symmetrize the square input matrix in, or a stack
% of images.
% If a second parameter angBW is given, it is the order of angular Fourier
% components that is kept.  (Default is 0, meaning complete circular
% averaging.)

[n n1 nim]=size(in);

nr=ceil(n/2);
nt=2*n;
ctr=ceil((n+1)/2);
f=ifftshift(-nt/2:nt/2-1);  % frequencies for filter

out=in;
for i=1:nim
    ps=ToPolar(in(:,:,i),nr,2*n,1,ctr,ctr);
%     ps=gridToPolar(in(:,:,i),2*n);
    if nargin<2
        sf=mean(ps,2);
    else
        H=abs(f)<=angBW;
        sf=real(ifft(fft(ps,nt,2).*repmat(H,nr,1),nt,2));
    end;
    out(:,:,i)=ToRect(sf,n,1,[ctr ctr]);
%     out(:,:,i)=gridToRect(sf);
end;
