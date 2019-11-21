% modify AE from the one with J ambiguity to the one without.
% 0 = AE * [Q(:); Y(:)]; from AE * [Q(:); Y(5 of elements)]

load('AE.mat');
AE_short = AE;
qsiz = 16;
ysiz = 9;
siz0 = 8;

AE_noJ = zeros(qsiz+1,qsiz+ysiz); % constraints with 16 Q elements and 1 trace
                        % variable with 16 Q elements and 9 Y elements
AE_noJ(1:(1+siz0),1:qsiz) = AE_short(1:(1+siz0),1:qsiz);

AE_noJ(1:(1+siz0),qsiz+5) = AE_short(1:(1+siz0),qsiz+1);
AE_noJ(1:(1+siz0),qsiz+1) = AE_short(1:(1+siz0),qsiz+2);
AE_noJ(1:(1+siz0),qsiz+3) = AE_short(1:(1+siz0),qsiz+3);
AE_noJ(1:(1+siz0),qsiz+7) = AE_short(1:(1+siz0),qsiz+4);
AE_noJ(1:(1+siz0),qsiz+9) = AE_short(1:(1+siz0),qsiz+5);

AE_noJ(1+siz0+1,2) = 4; % Q(2,1)
AE_noJ(1+siz0+1,qsiz+2) = 1; 
AE_noJ(1+siz0+1,qsiz+4) = -1;
AE_noJ(1+siz0+2,5) = 1; % Q(1,2)
AE_noJ(1+siz0+2,2) = -1; 
AE_noJ(1+siz0+3,12) = 4; % Q(4,3)
AE_noJ(1+siz0+3,qsiz+2) = -1; 
AE_noJ(1+siz0+3,qsiz+4) = -1;
AE_noJ(1+siz0+4,15) = 1; % Q(3,4)
AE_noJ(1+siz0+4,12) = -1;
AE_noJ(1+siz0+5,3) = 4; % Q(3,1)
AE_noJ(1+siz0+5,qsiz+8) = 1; 
AE_noJ(1+siz0+5,qsiz+6) = -1;
AE_noJ(1+siz0+6,9) = 1; % Q(1,3)
AE_noJ(1+siz0+6,3) = -1; 
AE_noJ(1+siz0+7,8) = 4; % Q(4,2)
AE_noJ(1+siz0+7,qsiz+6) = 1; 
AE_noJ(1+siz0+7,qsiz+8) = 1;
AE_noJ(1+siz0+8,14) = 1; % Q(2,4)
AE_noJ(1+siz0+8,8) = -1;

AEAEtInv_noJ = pinv(AE_noJ*AE_noJ');

save('AE_noJ.mat','AE_noJ','AEAEtInv_noJ');

