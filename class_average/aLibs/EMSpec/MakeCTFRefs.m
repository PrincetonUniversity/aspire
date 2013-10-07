function [refs refsz]=MakeCTFRefs(n,res,P,msk)
% function [refs refsz]=MakeCTFRefs(n,res,P,msk)
% Compute a huge array of CTF functions for correlation with a power
% spetrum, to implement the algorithm of CTFFIND3.
% The structure P now has either scalars or vectors for the fields
% P.defocus, P.theta, P.alpha and P.B.
% There is a new field that also can be a vector, P.deltadef which is the
% amplitude of the astigmatism.  The defocus values def1 and def2 are given
% by P.defocus +/- P.deltadef.
% refs is an array of size n/2 x n x prod(refsize) where the sequence of
% subscripts are given by refsz=[nB nalpha ndef ndeldef ntheta]

defvals=P.defocus;
ndef=numel(defvals);

ddefvals=P.deltadef;
nddef=numel(ddefvals);

thetavals=P.theta;
ntheta=numel(thetavals);

alphavals=P.alpha;
nalpha=numel(alphavals);

Bvals=P.B;
nB=numel(Bvals);

% Construct radius and theta with half of x values.
[X,Y]=ndgrid(0:n/2-1,-n/2:n/2-1);
r2=(X.^2+Y.^2);
theta=atan2(Y,X);
f0=1/(n*res);
refsz=[nB, nalpha, ndef, nddef, ntheta];
refs=single(zeros(n/2,n,prod(refsz)));

chi4=pi/2*P.Cs*P.lambda^3*1e7*f0^4*r2.^2; % f^4 term of chi
k2=-pi*P.lambda*1e4*f0^2;  % constant for f^2 term
ii=0;

% pre-computer the envelope function
env=zeros(n/2,n,nB);
for m=1:nB
    env(:,:,m)=exp(-r2*f0^2*Bvals(m));
end;

for k=1:ntheta
    for j=1:nddef;
        for i=1:ndef
            % compute defocus values in the plane.
            df=defvals(i)+ddefvals(j)*cos(2*(thetavals(k)-theta));
            for l=1:nalpha
                % compute full CTF except for envelope
                c0=sin(k2*df.*r2+chi4-alphavals(l));
                for m=1:nB
                    ii=ii+1;
                    % c=msk.*c0.*exp(-r2*f0^2*Bvals(m));
                    %   no, look up the envelope function
                    c=(msk.*c0.*env(:,:,m)).^2;
                    refs(:,:,ii)=c/sqrt(c(:)'*c(:));  % normalize variance.
                end;
            end;
        end;
    end;
end;
