%% convert from cartesian volume to spherical fourier coefficients
load('cleanrib.mat')
n = size(volref,1);
%n = n*sqrt(n);
fac = 2;
az = linspace(-pi,pi,n*fac);
el = linspace(-pi/2,pi/2,n);
[az1,el1] = ndgrid(az,el);
[x,y,z] = sph2cart(az1,el1,ones(n*fac,n));
omega = [x(:) y(:) z(:)];
rate = 16;
st1 = n/2-rate/2+1;
st2 = n/2+rate/2;
vol1 = volref(st1:st2,st1:st2,st1:st2);
%f1 = nufft_3d(volref(st1:st2,st1:st2,st1:st2),omega,'single',rate);
Nd = [rate rate rate];
Jd = Nd/2;
Kd = Nd*2;
n_shift = Nd/2;
st = nufft_init(omega,Nd,Jd,Kd,n_shift,'table', 2^11, 'minmax:kb');
f1 = nufft(vol1,st);

%f11 = reshape(f1,n*n,1);

%% convert to spherical harmonics representation
lat = el1(:)/pi*180;
lon = az1(:)/pi*180;
[lmcosiY,dwY] = xyz2plm(f1,24,'im',lat,lon);%Size 625 / 650
figure()
plot(abs(lmcosiY(:,3))); hold on;
plot(abs(lmcosiY(:,4))); 
legend('Real part of SH coefficients','complex part of SH coefficients');
xlabel('l,m index');
ylabel('magnitude');
title('Energy of SH coefficients, clean 16^3 volume');

[r,~,~,~] = plm2xyz(lmcosiY,lat,lon);
err = f1-r;
e = mean(abs(err).^2);
p = mean(abs(f1).^2);
display(['Ratio between avarege squared error and image power: ', num2str(e/p)])

[u1,s1,v1] = svd(reshape(f1,n,n));
[u2,s2,v2] = svd(reshape(r,n,n));
figure()
plot(diag(s1),'o-');
hold on
plot(diag(s2),'+-');
legend('Fourier coefficients','Inversed from spherical harmonics');
title('spectrum before and after Spherical Harmonics transform');

%% rotate the volume and align
vol2 = zeros(rate,rate,rate);
%vol1 = volref(st1:st2,st1:st2,st1:st2);

% rotate volume clockwise around y axis => rotate coordinate
% counterclockwise
for i = 1:rate
    for j = 1:rate
        for k = 1:rate
            vol2(i,j,k) = vol1(i,rate+1-k,j);
        end
    end
end
figure()
subplot(2,2,1); imagesc(abs(vol1(:,:,1))); colormap(gray); title('Vol1 frame1')
subplot(2,2,2); imagesc(abs(vol2(:,:,1))); colormap(gray); title('Vol2 frame1')
subplot(2,2,3); imagesc(abs(vol1(:,:,8))); colormap(gray); title('Vol1 frame8')
subplot(2,2,4); imagesc(abs(vol2(:,:,8))); colormap(gray); title('Vol2 frame8')


f2 = nufft(vol2,st);

[lmcosiY2,dwY2] = xyz2plm(f2,24,'im',lat,lon);
figure()
plot(abs(lmcosiY2(:,3))); hold on;
plot(abs(lmcosiY2(:,4))); 
legend('Real part of SH coefficients','complex part of SH coefficients');
xlabel('l,m index');
ylabel('magnitude');
title('Energy of SH coefficients, clean 16^3 volume');

[r2,~,~,Plm2] = plm2xyz(lmcosiY2,lat,lon);
% [u21,s21,v21] = svd(reshape(f2,n,n));
% [u22,s22,v22] = svd(reshape(r2,n,n));
% figure()
% plot(diag(s21),'o-');
% hold on
% plot(diag(s22),'+-');
% legend('Fourier coefficients','Inversed from spherical harmonics');
% title('spectrum before and after Spherical Harmonics transform, vol2');

%% rotate Fourier coefficient to check: optimization worked with rotating
%Fourier coefficient -> problem happened when rotating the volume &
%transform to Fourier -> checked again for NUFFT! with rot_est = [0 -1 0; 1
%0 0; 0 0 1]

az2 = linspace(-pi/2,pi/2*3,n);
el2 = linspace(-pi/2,pi/2,n);
[az2,el2] = ndgrid(az2,el2);
[x2,y2,z2] = sph2cart(az2,el2,ones(n,n));
omega2 = [x2(:) y2(:) z2(:)];

st22 = nufft_init(omega2,Nd,Jd,Kd,n_shift,'table', 2^11, 'minmax:kb');
f22 = nufft(vol1,st22);

lat2 = el2(:)/pi*180;
lon2 = az2(:)/pi*180;
[lmcosiY3,dwY3] = xyz2plm(f22,62,'im',lat,lon);
[r3,~,~,Plm] = plm2xyz(lmcosiY3,lat2,lon2);
[u3,s3,v3] = svd(reshape(f22,n,n));
[u32,s32,v32] = svd(reshape(r3,n,n));
figure()
plot(diag(s3),'o-');
hold on
plot(diag(s32),'+-');
legend('Fourier coefficients','Inversed from spherical harmonics');
title('spectrum before and after Spherical Harmonics transform');

%% align in SH coefficients
opt = zeros(1,3);
score = 10;
I = zeros(1,9*9*5);
A = cell(1,9*9*5);
iter = 0;
for alpha = -180:45:135
    for beta = -90:45:45
        for gamma = -135:45:180
            [lmcosiYr,~,~] = plm2rot(lmcosiY,alpha,beta,gamma);
            i = norm(lmcosiYr(:,3:4)-lmcosiY2(:,3:4));
            disp([alpha,beta,gamma,i]);
            iter = iter+1;
            I(iter) = i;
            A{1,iter} = [alpha beta gamma];
            if i < score
                score = i;
                opt = [alpha beta gamma];
            end
        end
    end
end
I = I(1:256);
m = find(I==min(I(:)));
A{1,I==min(I)}
rot_est = EA_MAT_SOFT(opt/180*pi);
figure()
plot(I)

            