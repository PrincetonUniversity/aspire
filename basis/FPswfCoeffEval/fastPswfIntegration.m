function [ coeffVecQuadFast ] = fastPswfIntegration( images_nufft, c, L, numAngularPts, ang_freq, radialQuadPts, quadRuleRadialWts, PSWF_radial_quad, realFlag)
% This function computes the PSWF expansion coefficients based on the NFFT
% of the images and the quadrature nodes and weights supplied.

nImages = size(images_nufft,2);
if (realFlag==1)
    r_quad_idx = cumsum([1;numAngularPts/2]);
else
    r_quad_idx = cumsum([1;numAngularPts]);
end
N_max = max(ang_freq)+1;
r_N_eval_mat = zeros(numel(radialQuadPts),N_max,nImages);
for i = 1:numel(radialQuadPts)
    if (realFlag==1)
        curr_r_mat = images_nufft(r_quad_idx(i):(r_quad_idx(i)+numAngularPts(i)/2-1),:);
        curr_r_mat = [curr_r_mat; conj(curr_r_mat)];
    else
        curr_r_mat = images_nufft(r_quad_idx(i):(r_quad_idx(i)+numAngularPts(i)-1),:);
    end
    angularEval = fft(curr_r_mat)*quadRuleRadialWts(i);
    if (numAngularPts(i)<N_max)
        angularEval = repmat(angularEval,ceil(N_max/numAngularPts(i)),1);
    end
    r_N_eval_mat(i,:,:) = angularEval(1:N_max,:);    
end

numelForN = zeros(1,N_max);
for i = 0:N_max-1
    numelForN(i+1) = nnz(ang_freq==i);
end
idxForN = cumsum([1,numelForN]);

r_N_eval_mat = reshape(r_N_eval_mat,numel(radialQuadPts)*N_max,nImages);
coeffVecQuadFast = zeros(numel(ang_freq),nImages);
for i = 0:N_max-1
    currBlkR = c/(2*pi*L)*PSWF_radial_quad(:,idxForN(i+1)-1+(1:numelForN(i+1))).'; 
    coeffVecQuadFast(idxForN(i+1)-1+(1:numelForN(i+1)),:) = currBlkR * r_N_eval_mat(i*numel(radialQuadPts)+(1:numel(radialQuadPts)),:);    
end


end

