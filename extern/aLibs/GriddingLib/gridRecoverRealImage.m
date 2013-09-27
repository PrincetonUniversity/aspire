function m=gridRecoverRealImage(P, postcomp)
% function m=gridRecoverRealImage(P, [postcomp]);
% Transform the padded FT back to a real image or volume.
% If postcompensation is desired, use gridMakePreComp to create the column
% vector postcomp.  Default is postcommp=1.
% The output is masked by fuzzydisc with a radius n/2-2.
% This function handles 1d, 2d or 3d inputs.

if nargin<2
    postcomp=1;
else
    postcomp=postcomp(:);  % force it to be a column vector
end;

dimension=ndims(P.PadFT);
sizes=size(P.PadFT);
if sizes(2)==1
    dimension=1;
end;

np=P.np;
sp=P.sp;
n=P.n;
% There are two possibilities: P.PadFT could be np or np1 in size.
% to handle this, we compute the shift based on the actual P.PadFT size.
sp1=(sizes(1)-np)/2;
maskr=n/2+1;  % mask radius

switch dimension
    case 3
        rft=FromCAS(P.PadFT(sp1+1:sp1+np, sp1+1:sp1+np, sp1+1:sp1+np));
        rm=fftshift(real(ifftn(fftshift(rft))));
        if numel(postcomp)>1
            pc2=kron(postcomp,postcomp');
            rm=rm.*reshape(kron(pc2(:),postcomp'),[np np np]);  % perform the 3D compensation
        end;
        rmcm=rm.*fuzzymask(np,3,maskr,3);  % mask out the outlying part
        m=rmcm(1+sp:n+sp, 1+sp:n+sp, 1+sp:n+sp);  % extract the un-padded result.
    case 2
        mask=fuzzymask(np,2,maskr,3);
        %     postcomp=postcomp'./max(.01,masksum);
        rft=FromCAS(P.PadFT(sp1+1:sp1+np,sp1+1:sp1+np));
        rm=fftshift(real(ifftn(fftshift(rft)))); % Get the np x np reconstruction.
        rm=rm.*kron(postcomp,postcomp');  % perform the post-compensation
        rmcm=rm.*mask;  % mask out the outlying part
        m=rmcm(1+sp:n+sp,1+sp:n+sp);  % extract the un-padded result.
    case 1
        rft=FromCAS(P.PadFT(sp1+1:sp1+np));
        rm=fftshift(real(ifft(fftshift(rft)))); % Get the np x np reconstruction.
        rm=rm.*postcomp;  % perform the post-compensation
        mask=fuzzymask(np,1,maskr,3);
        rmcm=rm.*mask;  % mask out the outlying part
        m=rmcm(1+sp:n+sp);  % extract the un-padded result.
end;
