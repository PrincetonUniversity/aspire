% Compare execution speeds of the old PFT code (cryo_pft_legacy) and its
% implmenetation based on CIMS NUFFT (cryo_pft).
% 
% Yoel Shkolnisky, April 2017.

n_vals=32:16:512;
n_vals=n_vals+1; % For odd sizes
timing=zeros(numel(n_vals),2);
for i=1:numel(n_vals)
    n=n_vals(i);
    fprintf('Testing n=%d',n);
    n_projs=100;
    A=rand(n,n,n_projs);
    n_r=ceil(n/2);
    L=360;
    t1=tic;pf1=cryo_pft_legacy(A,n_r,L,'single');t1=toc(t1);
    t2=tic;pf2=cryo_pft(A,n_r,L,'single');t2=toc(t2);
    err=norm(pf1(:)-pf2(:))/norm(pf1(:));
    fprintf( '   err=%5.3e  t(legacy)/t(cims)=%5.3f\n',err,t1/t2);
    if  err>1.0e-5
        error('error too large');
    end
    timing(i,1)=t1; timing(i,2)=t2;
end
plot(n_vals,timing(:,1)./timing(:,2),'-o')
hold on;
plot(n_vals,ones(numel(n_vals),1),'r','LineWidth',2); % no gain line
hold off;

%disp([n_vals.' timing(:,1)./timing(:,2)]);