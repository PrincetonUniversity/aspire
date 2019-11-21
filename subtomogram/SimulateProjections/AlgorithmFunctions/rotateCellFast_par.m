function [ CR ] = rotateCellFast( cell1, angles, Wp, Wm )
% fast rotation of angular radial representation of a volume

% cell1: cell1{ll+1,1} represents the l-th order coefficient where m stands
% for column index, and radial index n stands for row index
% angles: 3*N vector representing Euler ZYZ angles of the desired rotations

% N = 10000;
% angles = rand(3,N)*2*pi;
% angles(2,:) = acos(rand(1,N)*2-1);
maxL = size(cell1,1)-1;
N = size(angles,2);

CR = cell(maxL+1,1);

parfor ll = 0:maxL
    atemp = reshape(transpose(cell1{ll+1,1}),2*ll+1,1,[]);
    atemp = bsxfun(@times, exp(-1i*(-ll:ll)'*angles(1,:)), atemp);
    atemp = reshape(atemp, 2*ll+1,[]);
    atemp = Wp{1,ll+1}*atemp;
    atemp = reshape(atemp, 2*ll+1, N, []);
    atemp = bsxfun(@times, exp(-1i*(-ll:ll)'*angles(2,:)), atemp);
    atemp = reshape(atemp, 2*ll+1, []);
    atemp = Wm{1,ll+1}*atemp;
    atemp = reshape(atemp, 2*ll+1, N, []);
    atemp = bsxfun(@times, exp(-1i*(-ll:ll)'*angles(3,:)), atemp);
    CR{ll+1,1} = permute(atemp(ll+1:end,:,:), [3,1,2]);
end

end

% obtain N random rotations for the volume cell1
% N = 2;
% angles = rand(3,N)*2*pi;
% angles(2,:) = acos(rand(1,N)*2-1);
% maxL = size(cell1,1)-1;
% 
% angles = [0,0,0; pi/2, pi/4, pi/6]';%angles(:,1)';
% 
% AN1 = cell(maxL+1,1);
% tic
% for ll = 0:maxL
%     Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
%     %[T,Tinv] = realY_to_complexY(ll);
%     %Wl = real(Tinv*Wl*T);
%     AN1{ll+1,1} = Wl'*transpose(cell1{ll+1});%cell1{ll+1}*conj(Wl); % give complex numbers, use real Y?
% end
% t_orig = toc
% 
% tic
% Ctemp = cell(1,maxL+1);
% Ctemp(:) = {[-pi/2, pi/2, pi/2]}; % Y axis -> Z axis
% Wp = cellfun(@wignerd, num2cell(0:maxL), Ctemp, 'UniformOutput', false);
% Ctemp(:) = {[-pi/2, -pi/2, pi/2]}; % Z' axis -> Y' axis
% Wm = cellfun(@wignerd, num2cell(0:maxL), Ctemp, 'UniformOutput', false);
% t_2pi = toc
% 
% % tic
% % for ll = 0:maxL % checked to eps: Ry(beta) = Rx(-pi/2) Rz(beta) Rx (pi/2)
% %     D{ll+1,1} = norm(wignerd(ll,[angle(1) 0 0])*...
% %         Wm{1,ll+1}*wignerd(ll,[angle(2) 0 0])*Wp{1,ll+1}...
% %         *wignerd(ll,[angle(3) 0 0]) - wignerd(ll,angle)); % from ll to -ll, default uses the - in exponent
% % end
% % t_diff = toc
% 
% W = cell(maxL+1,1);
% tic
% for ll = 0:maxL % checked to eps: Ry(beta) = Rx(-pi/2) Rz(beta) Rx (pi/2)
%     W{ll+1,1} = diag(exp(-1i*(-ll:ll)*angle(1)))*...
%         conj(Wm{1,ll+1})*diag(exp(-1i*(-ll:ll)*angle(2)))*conj(Wp{1,ll+1})...
%         *diag(exp(-1i*(-ll:ll)*angle(3))); % from ll to -ll, default uses the - in exponent
% end
% t_decomp = toc
% 
% W2 = cell(maxL+1,1);
% tic
% for ll = 0:maxL % checked to eps: Ry(beta) = Rx(-pi/2) Rz(beta) Rx (pi/2)
%     W2{ll+1,1} = diag(exp(-1i*(-ll:ll)*angle(3)))*...
%         conj(Wp{1,ll+1})*diag(exp(-1i*(-ll:ll)*angle(2)))*conj(Wm{1,ll+1})...
%         *diag(exp(-1i*(-ll:ll)*angle(1))); % from ll to -ll, default uses the - in exponent
% end
% t_decomp = toc
% 
% tic
% AN2 = cellfun(@(x,y) y*transpose(x), cell1, W2, 'UniformOutput', false);
% t_mult = toc;
% 
% D = cellfun(@(x,y) norm(x-y), AN1, AN2, 'UniformOutput', false);
% 
% CR = cell(maxL+1,1);
% tic
% for ll = 0:maxL
%     atemp = reshape(transpose(cell1{ll+1,1}),2*ll+1,1,[]);
%     atemp = bsxfun(@times, exp(-1i*(-ll:ll)'*angles(1,:)), atemp);%exp(-1i*(-ll:ll)'*angles(1,:))
%     atemp = reshape(atemp, 2*ll+1,[]);
%     atemp = conj(Wm{1,ll+1})*atemp;
%     atemp = reshape(atemp, 2*ll+1, N, []);
%     atemp = bsxfun(@times, exp(-1i*(-ll:ll)'*angles(2,:)), atemp);
%     atemp = reshape(atemp, 2*ll+1, []);
%     atemp = conj(Wp{1,ll+1})*atemp;
%     atemp = reshape(atemp, 2*ll+1, N, []);
%     atemp = bsxfun(@times, exp(-1i*(-ll:ll)'*angles(3,:)), atemp);
%     CR{ll+1,1} = permute(atemp, [1,3,2]);
%     %CR{ll+1,1} = reshape(atemp, 2*ll+1, []);
% end
% toc
% 
% D2 = cellfun(@(x,y) norm(x-y), AN1, CR, 'UniformOutput', false);
