cd /u/liuyuan/Documents/MATLAB/Subtomo_average/aspire
initpath;
addpath(genpath('/u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE/subtomogram_yuan'));
addpath([getenv('HOME') '/local/share/nfft/matlab/nfft']);

N=5;
Eavg = zeros(3,2,5);
Tavg = zeros(3,3,5);

for j = 4:8
    siz = 2^j;
    E = zeros(3,2,N);
    T = zeros(3,3,N);
    nj=siz^3;
    n=siz;
    for i = 1:N
        [e_iter,t_iter] = nufft_compare(nj,n);
        E(:,:,i) = e_iter;
        T(:,:,i) = t_iter;
    end
    Tavg(:,:,j-3) = mean(T,3);
    Eavg(:,:,j-3) = mean(E,3);
    save('nufft_compare4.mat')
end

%%
% load('nufft_compare.mat');
% n = 2.^(4:6);
% m = n.^3;
% comp_g = n.^3.*log(n)+30^3*m;
% comp_k = n.^3.*log(n)+6^3*m;
% comp_y = (2*n).^3.*log(n)+log(10^16)*m;
% t = zeros(1,3);
% figure(1)
% subplot(3,1,1);
% plot(comp_g,reshape(Tavg(1,1,1:3),1,3),'o-',comp_g,reshape(Tavg(1,2,1:3),1,3),'o-');
% xlabel('Order of algorithm complexity');
% ylabel('Computation time (s)');
% title('Algorithm complexity for Greengard 3D NUFFT with 16^3, 32^3, 64^3 volumes',...
%     'FontSize',13);
% legend('Forward transform','Inverse transform','Location','southeast');
% subplot(3,1,2);
% plot(comp_k,reshape(Tavg(2,1,1:3),1,3),'o-',comp_k,reshape(Tavg(2,2,1:3),1,3),'o-');%6^3*m,reshape(Tavg(2,3,1:3),1,3));
% xlabel('Order of algorithm complexity');
% ylabel('Computation time (s)');
% title('Algorithm complexity for Keiner 3D NFFT with 16^3, 32^3, 64^3 volumes',...
%     'FontSize',13);
% legend('Forward transform','Inverse transform','Location','southeast');
% subplot(3,1,3);
% plot(comp_k,reshape(Tavg(3,1,1:3),1,3),'o-',comp_k,reshape(Tavg(3,2,1:3),1,3),'o-');%6^3*m,reshape(Tavg(2,3,1:3),1,3));
% xlabel('Order of algorithm complexity');
% ylabel('Computation time (s)');
% title('Algorithm complexity for Yoel 3D NUFFT with 16^3, 32^3, 64^3 volumes',...
%     'FontSize',13);
% legend('Forward transform','Inverse transform','Location','east');
% 

