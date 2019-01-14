
%Input:  P = Fourier transform of projection image
%        B = Band limit for speed and accuracy
%Output: m = A measure of how much P is an equator
%        s = angle (in degrees) of in_plane rotation
%        scl = The values on the approximated self common line of the image
function [m,scl,s_out]=cryo_eqMeasure(P,B)

%P is nrxL2
L2=size(P,2);
L=floor(L2/2);
ninety_deg=floor(L2/4);
m=inf;
scl=P(:,size(P,2));
for s=0:L-1
    %Rotate image by angle s
    Ps=cat(2,P(1:B,s+1:end),P(1:B,1:s));
    %Test equator property
    tmp=abs(Ps(2:B,1:L/2)-flip(conj(Ps(2:B,L/2+2:L+1)),2))./abs(Ps(2:B,L/2+1:L));
    tmp=(mean(mean(tmp,1)));%+...
        %mean(abs((Ps(2:B,L/2+1)-Ps(2:B,L/2+1+L)))./abs(Ps(2:B,L/2+1+L))))/2;
    if tmp<m
        m=tmp;
        scl=Ps(:,L/2+1+s+ninety_deg);
        s_out=s;
    end
end
