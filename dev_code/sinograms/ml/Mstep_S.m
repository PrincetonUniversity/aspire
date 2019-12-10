function S=Mstep_S(X,Z,Zi,S,ht,K,KNN,shifts,max_shift,shift_step) 

rmax=size(X,1);
rk2=(0:(rmax-1))';
n_shifts=ceil(2*max_shift/shift_step+1); % Number of shifts to try.
shift_phases=zeros(rmax,n_shifts);
for shiftidx=1:n_shifts
    shift=-max_shift+(shiftidx-1)*shift_step;
    shift_phases(:,shiftidx)=exp(2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
end

% for k=1:K
%     ind=find(Zi==k);
%     if isempty(ind)
%         continue;
%     end
%     nl=ceil(ind/KNN);
%     Unlk=shift_phases(:,shifts(ind));  
%     Znlk=Z(ind);
%     S(:,k)=sum(bsxfun(@times,(ht(nl).*Znlk)',Unlk.*X(:,nl)),2)/sum(ht(nl).^2.*Znlk);
% end
[Zis,ZiI]=sort(Zi(:));
% tic;
% for k=1:K
%     [lower_index,upper_index] = bsearch(Zis,k,k);
%     if isempty(lower_index) || isempty(upper_index)
%         continue;
%     end
%     ind=ZiI(lower_index:upper_index);
%         nl=ceil(ind/KNN);
%     Unlk=shift_phases(:,shifts(ind));  
%     Znlk=Z(ind);
%     
%     % XXX Why is ht of length N*L and not of length N?
%     % XXX It seems that you access
%     S(:,k)=sum(bsxfun(@times,(ht(nl).*Znlk)',Unlk.*X(:,nl)),2)/sum(ht(nl).^2.*Znlk);
% end
% toc

%S2=zeros(size(S));
%tic;

reset(gpuDevice); % To clear memory. 
ght=gpuArray(single(ht));
gZ=gpuArray(single(Z));
gdenom=bsxfun(@times,(ght.^2).',gZ);
gnom1=bsxfun(@times,ght.',gZ);
gnom1=conj(gnom1);
%gshift_phases=gpuArray(single(shift_phases));
gS=gpuArray(single(zeros(size(S))));
for k=1:K
    [lower_index,upper_index] = bsearch(Zis,k,k);
    if isempty(lower_index) || isempty(upper_index)
        continue;
    end
    ind=ZiI(lower_index:upper_index);
    nl=ceil(ind/KNN);
    %gUnlk=gshift_phases(:,shifts(ind));  
    %Znlk=Z(ind);
    
    % XXX Why is ht of length N*L and not of length N?
    % XXX Why do you have imagninary valyes in the denominator?
    tmp1=shift_phases(:,shifts(ind)).*X(:,nl);
    gtmp1=gpuArray(tmp1);
    %gtmp1=gUnlk.*X(:,nl);
    gS(:,k)=sum(bsxfun(@times,(gnom1(ind)).',gtmp1),2)/sum(gdenom(ind));
end
S=gather(gS);
assert(~any(isnan(S(:))))    % S must not contain NaNs.
%toc
%aaa=1;

