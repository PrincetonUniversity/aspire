function OUTPUT=fastrotate_gpu(INPUT,phi,M,precision)
%FASTROTATE_GPU Rotate a 2D image.
%
% GPU implementation of fastrotate.
% See fastrotate for information.
%
%Yoel Shkolnisky, April 2014.

if ~exist('precision','var')
    precision='single';
end

[SzX SzY ~] =size(INPUT);

comp=0; % Determine if we need to compute the interpolation tables
if ~exist('M','var')
    comp=1;
elseif ~isstruct(M)
    comp=1;
elseif ~isfield(M,'Mx') || ~isfield(M,'My');
    comp=1;
else
    Mx=M.Mx;
    My=M.My;
end

if comp
    % Precompte My and Mx
    M=fastrotateprecomp(SzX,SzY,phi);
    Mx=M.Mx;
    My=M.My;
end

if strcmpi('precision','single')
    Mx=single(Mx);
    My=single(My);
    INPUT=single(INPUT);
end

gMx=gpuArray(Mx);
gMy=gpuArray(My);
gINPUT=gpuArray(INPUT); % Set precision
gOUTPUT=gpuArray.zeros(size(INPUT),precision);

    

parfor k=1:size(INPUT,3)
    
    % Rotate by multiples of 90 degrees.
    switch M.mult90
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
OUTPUT=double(gather(gOUTPUT));
