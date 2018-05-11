function gOUTPUT=fastrotate3daux_gpu(gINPUT,gMx,gMy,mult90,precision)
%
% Auxiliry function for fastrotate3d_gpu
%
% Yoel Shkolnisky, April 2014.

gOUTPUT=gpuArray.zeros(size(gINPUT),precision);

for k=1:size(gINPUT,3)
    
    % Rotate by multiples of 90 degrees.
    switch mult90
        case 0
        case 1
            gINPUT(:,:,k)=rot90_fastrotate(gINPUT(:,:,k));
        case 2
            gINPUT(:,:,k)=rot180(gINPUT(:,:,k));
        case 3
            gINPUT(:,:,k)=rot270(gINPUT(:,:,k));
        otherwise
            error('Invalid value for mult90');
    end
    
    %FIRST PASS
    % for x=1:SzX,
    %     spinput=fft(INPUT(x,:));
    %     spinput=spinput.*(My(:,x)).';
    %     OUTPUT(x,:)=real(ifft(spinput));
    % end
    
    gspinput=fft(gINPUT(:,:,k),[],2);
    gspinput=gspinput.*gMy;
    gOUTPUT(:,:,k)=real(ifft(gspinput,[],2));
    
    % SECOND PASS
    % for y=1:SzY,
    %     spinput=fft(OUTPUT(:,y));
    %     spinput=spinput.*Mx(:,y);
    %     OUTPUT(:,y)=real(ifft(spinput));
    % end
    
    gspinput=fft(gOUTPUT(:,:,k));
    gspinput=gspinput.*gMx;
    gOUTPUT(:,:,k)=real(ifft(gspinput));
    
    %THIRD PASS
    % for x=1:SzX,
    %     spinput=fft(OUTPUT(x,:));
    %     spinput=spinput.*(My(:,x)).';
    %     OUTPUT(x,:)=real(ifft(spinput));
    % end
    
    gspinput=fft(gOUTPUT(:,:,k),[],2);
    gspinput=gspinput.*gMy;
    gOUTPUT(:,:,k)=real(ifft(gspinput,[],2));
    
end

