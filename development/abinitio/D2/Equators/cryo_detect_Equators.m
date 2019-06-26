
%Input:  Pf = Fourier transformed images
%        cutoff = percentile for cuttoff
%        B = Band limit for cryo_eqMeasure
%Output: 
function [eq_idx,eqm,eq_class,rm]=cryo_detect_Equators(pf,B,cutoff)

nImages=size(pf,3);
%nr=size(pf,1);
scls=zeros(B,nImages);
eqm=zeros(nImages,1);
for i=1:nImages
    P=pf(:,:,i);
    [eqm(i),scls(:,i),s]=cryo_eqMeasure(P,B);
end

histogram(eqm,10);
eq_idx=eqm<=prctile(eqm,cutoff);%find(eqm<=prctile(eqm,cutoff));
scls=scls(2:end,eq_idx);
[eq_class,C,~,D]=kmeans(real(scls'),3);
rm=log10(abs(real(scls)./imag(scls)));




