function [ AN1 ] = rotateCell( A, angle )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

maxL = size(A,1)-1;
AN1 = cell(maxL+1,1);
for ll = 0:maxL
    Wl = wignerd(ll,angle); % from ll to -ll, default uses the - in exponent
    %[T,Tinv] = realY_to_complexY(ll);
    %Wl = real(Tinv*Wl*T);
    al = A{ll+1}*conj(Wl); % give complex numbers, use real Y?
    AN1{ll+1,1} = al;
end

end

%%


%cellR = rotateCell(cell1,[pi/3,pi/4,pi/5]);
angle = [pi/3,pi/4,pi/5];
maxL = size(cell1,1)-1;

tic
Ctemp = cell(1,maxL+1);
Ctemp(:) = {[-pi/2, pi/2, pi/2]}; % Y axis -> Z axis
Wp = cellfun(@wignerd, num2cell(0:maxL), Ctemp, 'UniformOutput', false);
Ctemp(:) = {[-pi/2, -pi/2, pi/2]}; % Z' axis -> Y' axis
Wm = cellfun(@wignerd, num2cell(0:maxL), Ctemp, 'UniformOutput', false);
toc

tic
for ll = 0:maxL % checked to eps: Ry(beta) = Rx(-pi/2) Rz(beta) Rx (pi/2)
    D{ll+1,1} = norm(wignerd(ll,[angle(1) 0 0])*...
        Wm{1,ll+1}*wignerd(ll,[angle(2) 0 0])*Wp{1,ll+1}...
        *wignerd(ll,[angle(3) 0 0]) - wignerd(ll,angle)); % from ll to -ll, default uses the - in exponent
end
toc
tic
AN1 = cellfun(@(x,y) x*conj(y), cell1, W, 'UniformOutput', false);
toc;