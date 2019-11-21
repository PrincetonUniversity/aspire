function [ Fnm, err ] = sht_cleanrib( f1, L, az1, el1 )
% Transform of the fourier coefficients to spherical harmonics expansion,
% with largest order L, latitude, and longitude specified.
% Output the positive part of the expansion coefficients, and the err when
% transforming back.

lat = el1(:)/pi*180;
lon = az1(:)/pi*180;

[lmcosiY,~] = xyz2plm(real(f1),L,'im',lat,lon);%Size 625 / 650
[lmcosiYc,~] = xyz2plm((f1-real(f1))*(-1i),L,'im',lat,lon);
% figure()
% plot(abs(lmcosiY(:,3))); hold on;
% plot(abs(lmcosiY(:,4))); 
% legend('Real part of SH coefficients','complex part of SH coefficients');
% xlabel('l,m index');
% ylabel('magnitude');
% title('Energy of SH coefficients, clean 32^3 volume, L58,k2');
[r,~,~,~] = plm2xyz(lmcosiY,lat,lon,L);
[rc,~,~,~] = plm2xyz(lmcosiYc,lat,lon,L);

err = norm(f1-r-rc*1i)/(norm(abs(f1))+norm(abs(r+rc*1i)))*2;

fnm = lmcosiY(:,3:4) + lmcosiYc(:,3:4);
Fnm = fnm(:,1)+fnm(:,2)*1i;

end

