function W = VesicleFromModel2(n,a,model,org)
%  function W = VesicleFromModel2(n,a,model,org)
% Version 2 uses interpolation of an oversampled 1d function for higher
% precision of high-frequency components.
% Given a density cross-section vector 'model', construct the density of a
% spherical vesicle in an image of size n, with radius a and center org.
% Let nm be numel(model). Then the density of the center of
% the membrane is taken to be model(nm/2+1) when nm is even, or
% model((nm+1)/2) when nm is odd; this is the same convention as fftshift.
% VesicleFromModel(n,a,ones(d,1),org) gives the same result as
% VesicleDensity(n,a,d,org).  n can be a 2-element vector to make a
% rectangular result.
% n=256; a=50; model=[0 1 0];
% org=[129 129];
% figure(1); SetGrayscale;
narrowShells=0;

tol = 1e-3;  % minimum density different to invoke calculation of a shell
ovs=8;  % oversampling factor

if numel(n)<2  % scalar means a square image
    n=[n n];
end;
ctr=ceil((n+1)/2);  % correct also for odd n
if nargin<4
    org=ctr;
end;

org=org(:)';  % make it a row vector
n1=[1 1]*2*ceil(a+numel(model)/2)+2;  % minimim square that bounds the vesicle.
if n1<min(n)  % We'll model the vesicle in a smaller square area
    ctr1=ceil((n1+1)/2);
    fracShift=org-round(org);
    org1=fracShift+ctr1;
else
    org1=org;
    n1=n;
end;

%% Create the band-limited, oversampled 1D function

modx=[0;model(:);0];  % pad with zeros
nm=numel(modx);
maxDensity=max(abs(modx));
mctr=ceil((nm+1)/2);

% Define the 1d radius
nr=ceil((a+mctr)*ovs)+1;
r1=(0:nr-1)/ovs;  % radius is taken in small steps
w1=zeros(1,nr);

rs=r1/ovs;
rp=r1+0.5/ovs;
rm=r1-0.5/ovs;

sp=rp.^2;
sm=rm.^2;
if narrowShells
    for i=2:nm-1
        %     We integrate shells 1/ovs pixel thick
        r0=i+a-mctr-.5/ovs;
        s0=r0^2;
        r1=r0+1/ovs;
        s1=r1^2;
        b=modx(i);
        % density of spherical shell
        w1=w1-ovs*b*real(rp.*sqrt(s0-sp)+s0*atan(rp./(sqrt(s0-sp)))...
                    -rm.*sqrt(s0-sm)-s0*atan(rm./(sqrt(s0-sm)))...
                    -rp.*sqrt(s1-sp)-s1*atan(rp./(sqrt(s1-sp)))...
                    +rm.*sqrt(s1-sm)+s1*atan(rm./(sqrt(s1-sm))));
    end;
else
    
    for i=2:nm
        %     We integrate shells one pixel thick
        r0=i+a-mctr-.5;
        s0=r0^2;
        b=modx(i-1)-modx(i);
        if r0>0 && abs(b)>tol*maxDensity
            % density of spherical shell
            w1=w1+b*real(rp.*sqrt(s0-sp)+s0*atan(rp./(sqrt(s0-sp)))...
                -rm.*sqrt(s0-sm)-s0*atan(rm./(sqrt(s0-sm))));
        end;
    end;
end;

% Filter and downsample the 1D profile
w2=[fliplr(w1(2:nr)) w1];  % symmetrize
nw=numel(w2);
nx=NextNiceNumber(numel(w2)*sqrt(2));
wx=ifftshift(Crop(w2,nx));  % pad the function and shift the origin
wxf=SharpFilt(wx,0.45/ovs,0.1/ovs);

%% Interpolate into the 2D field

r2=Radius(n1,org1)*ovs+1;
wi=ovs*interp1(wxf,r2(:),'spline');
wi=reshape(wi,n1);

if any(n1<n)
    W=ExtractImage(wi,round(org),n,1);  % insert the image into the larger one.
else
    W=wi;
end;

return

subplot(223);
plot([wx(40*ovs:55*ovs) wxf(40*ovs:55*ovs)]);


subplot(222);
% W=W;
imacs(W);
subplot(224);

W0 = VesicleFromModel(n,a,model,org);
q1=sect(W); q2=sect(W0);
plot([q1(75:85) q2(75:85)],'.-','markersize',10);

subplot(221);
semilogy([RadialPowerSpectrum(W) RadialPowerSpectrum(W0)]);
