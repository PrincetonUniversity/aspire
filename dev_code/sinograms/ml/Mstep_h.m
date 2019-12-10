function [h]=Mstep_h(X,Z,Zi,St,h,C,N,L)

% tic;
% for n=1:N*L
%     Sk=St(:,Zi(:,n));
%     Znlk=Z(:,n);
%     Xnl=X(:,n);
%     % XXX Why the real is only on the nominator? Why real at all?
%     h(n)=real((C*Xnl)'*Sk)*Znlk/(Znlk'*diag(Sk'*(C*Sk)));
% end
% toc

gC=gpuArray(single(C));
gSt=gpuArray(single(St));
gX=gpuArray(single(X));
gZ=gpuArray(single(Z));

gh=gpuArray(single(zeros(size(h))));
gCSt=gC*gSt;
gdn1=sum(conj(gSt).*gCSt); gdn1=gdn1.';
gCX=gC*gX; gCX=gCX';
for n=1:N*L
    gZnlk=gZ(:,n);
    gdenom=sum(conj(gZnlk).*gdn1(Zi(:,n),:));
    gh(n)=real(gCX(n,:)*gSt(:,Zi(:,n))*gZnlk/gdenom);
end
h=gather(gh);
assert(~any(isnan(h(:))))    % h must not contain NaNs.
%assert(norm(h(:)-h2(:))/norm(h(:))<1.0e-13);
