function ravg2d=cryo_radial_powerspect_2d(projs,pad,iswindow)

if nargin<2
    pad=0;
end

if nargin<3
    iswindow=0;
end

if pad==0
    if iswindow
        w=bartlett(size(projs,1));
        W=w*w.';
    else
        W=1;
    end
    projshat=zeros(size(projs));
    for k=1:size(projs,3)
        p=projs(:,:,k);
        p=cfft2(p).*W;
        p=(abs(p)).^2;
        projshat(:,:,k)=p;
    end
else
    n=size(projs,1);
    K=size(projs,3);
    
    if iswindow
        w=bartlett(2*n-1);
        W=w*w.';
    else
        W=1;
    end
    
    projshat=zeros(2*n-1,2*n-1,K);    
    for k=1:K
        p=zeros(2*n-1,2*n-1);
        p(n:2*n-1,n:2*n-1)=projs(:,:,k);
        p=cfft2(p).*W;
        p=(abs(p)).^2;
        projshat(:,:,k)=p;
    end
end
ravg2d=cryo_radial_average2d(mean(projshat,3));
orig2d=ceil((size(ravg2d+1)/2)); % Center of the 2D rdially averaged PSD of the projections.
ravg2d=ravg2d(orig2d(1):end,orig2d(2)); % Central ray of the PSD.
ravg2d=ravg2d./sum(ravg2d); % Normalize ray to mass 1.
