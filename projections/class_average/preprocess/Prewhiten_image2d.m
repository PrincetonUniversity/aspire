function [ proj ] = Prewhiten_image2d(proj, noise_response )
%Prewhiten Images with the power spectrum of the noise
%   noise_response is a 2d image
n=size(proj, 3);
L=size(proj, 1); 
l=floor(L/2);
K=size(noise_response, 1);
k=ceil(K/2);

%noise_response=noise_response(k-l:k+l, k-l:k+l);
%noise_response=noise_response/norm(noise_response(:));

filter=sqrt(noise_response);
filter=filter/norm(filter(:));

assert(norm(imag(filter(:)))<1.0e-14);
filter=real(filter);

assert(norm(filter-flipud(filter))<1.0e-13);
assert(norm(filter-fliplr(filter))<1.0e-13);

filter=(filter+flipud(filter))./2;
filter=(filter+fliplr(filter))./2;

nzidx=find(filter>1.0e-12);

%rmask=fuzzymask(K,2,0.45*K,0.05*K);

parfor i=1:n
    pp=zeros(K);
    pp(k-l:k+l, k-l:k+l)=proj(:, :, i);
    fp=cfft2(pp);
    p=zeros(size(fp));
    p(nzidx) = fp(nzidx)./filter(nzidx); 
    %p=p.*rmask;
    p2 = icfft2(p);
    assert(norm(imag(p2(:)))/norm(p2(:))<1.0e-13);
    p2 = p2(k-l:k+l, k-l:k+l);   
    proj(:, :, i)=real(p2);
end;

end

