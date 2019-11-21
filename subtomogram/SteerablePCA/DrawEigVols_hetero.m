% align the principal components, invert them back to the real domain,
% and draw central slices

% EigVols with 2^15 samples
load('eigVols.mat');
load('rib_basics.mat');
i = 1;
vecSt = U(:,i);
vecS = Us(:,i);
eul1 = alignComp(vecSt,vecS,siz_n);
% eul1(2) = -(180-eul1(2));
cellSR = rotateCell(vec2cell(vecSt,siz_n),[pi/4,0,0]);
eul1 = alignComp(cell2vec(cellSR),vecSt, siz_n);

%%
L = 16;
N = 2*L-1;
alpha = pi*(0:N)/L;
beta = pi*((0:N)*2+1)/(4*L);
[alpha,beta,gamma] = ndgrid(alpha,beta,alpha);
angles = [alpha(:) beta(:) gamma(:)];
A = vec2cell(vecSt,siz_n);
AN = cell(maxL+1, 1);
ANvec = zeros(size(Avec,1),(2*L)^3);
for i = 1:(2*L)^3
    angle = angles(i,:)';
    for ll = 0:maxL
        Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
        %[T,Tinv] = realY_to_complexY(ll);
        %Wl = real(Tinv*Wl*T);
        al = A{ll+1,1}*conj(Wl); % give complex numbers, use real Y?
        AN{ll+1,1} = al;
    end
    ANvec(:,i) = cell2vec(AN);
end
diff = zeros((2*L)^3,1);
for i = 1:(2*L)^3
    diff(i) = norm(ANvec(:,i)-vecS);
end
[~,ind] = min(diff);


%%
num = 5;
L = 9;
figure()
for i = 1:num
    vecSt = U(:,i);
    vecS = Us(:,i);
    eul2 = alignCompDirect(vecSt, vecS, siz_n, 9);
    cellSR = rotateCell(vec2cell(vecSt,siz_n),eul2);
    coeffSt = cell2coeff(cellSR);
    
    [ vol_hatSt, err, t_ift, ball ] = IFT_SB(R, c, L, coeffSt, F3d);
    subplot(2,num,i)
    threshF = quantile(vol_hatSt(:),0.1);
    Fv = isosurface(X,Y,Z,vol_hatSt,threshF);
    p = patch(Fv);
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1 1 1])
    view(3); 
    axis tight
    camlight 
    lighting gouraud
    title(['Steer PCA eigenvol,' num2str(i) num2str(sqrt(S(i)),'%10.3e')]);
end

for i = 1:num
    coeffS = cell2coeff(vec2cell(Us(:,i),siz_n));
    [ vol_hat1, err, t_ift, ball ] = IFT_SB(R, c, L, coeffS, F3d);
    subplot(2,num,num+i)
    threshF = quantile(vol_hat1(:),0.1);
    Fv = isosurface(X,Y,Z,vol_hat1,threshF);
    p = patch(Fv);
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1 1 1])
    view(3); 
    axis tight
    camlight 
    lighting gouraud
    title(['PCA eigenvol,' num2str(i) num2str(sqrt(Ss(i)),'%10.3e')]);
end

%%
num = 5;
mid = R+1;
figure()
for i = 1:num
    vecSt = U(:,i);
    vecS = Us(:,i);
    eul2 = alignCompDirect(vecSt, vecS, siz_n, 16);
    cellSR = rotateCell(vec2cell(vecSt,siz_n),eul2);
    coeffSt = cell2coeff(cellSR);
    
    [ vol_hatSt, err, t_ift, ball ] = IFT_SB(R, c, L, coeffSt, F3d);
    subplot(2,num,i)
    imagesc(vol_hatSt(:,:,mid));
    title(['Steer PCA eigenvol,' num2str(i) num2str(sqrt(S(i)),'%10.3d')]);
end

for i = 1:num
    cellS = vec2cell(Us(:,i), siz_n);
    coeffS = cell2coeff(cellS);
    [ vol_hat1, err, t_ift, ball ] = IFT_SB(R, c, L, coeffS, F3d);
    subplot(2,num,num+i)
    imagesc(vol_hat1(:,:,mid));
    title(['PCA eigenvol,' num2str(i) num2str(sqrt(Ss(i)),'%10.3d')]);
end
