function [proj_CTF, CTF, index]=  add_CTF_env_v6(proj, ndef, def1,def2,B, lambda, flag_CTF)

% Add CTF to a stack of images proj
% INPUTS:
% proj: stack of projections
% ndef: Number of defocus groups
% def1, def2: Min and max defocus
% B: envelope parameter
% lambda: Electron beam wavelength
% flag_CTF: 0 for no CTF, 1 for adding CTF
%
% OUTPUTS:
% proj_CTF: stack of projections with CTF applied
% CTF: stack of CTFs
% index: Index for the CTF for each image in the stack
% Tejal Bhamre: updated Oct 2015




dbstop if error

L=size(proj,1);
defocus=linspace(def1, def2, ndef);

    
if (flag_CTF==0)
       h=ones(L,L,ndef); 
       
elseif (flag_CTF==1)
       for i=1:ndef
       h(:, :, i) = CTF_old(L, 2.82, lambda, defocus(i), 2, B, 0.07);% generating CTF function using old CTF from Jane, new one compatible
       %with 80s data has more parameters and different units
       end  
elseif (nargin<7)
   error('Insufficient parameters') 
end

%% Add CTF 

% Generate projections with CTF
proj_CTF=zeros(size(proj));
index=round((1:size(proj,3))/size(proj,3)*ndef+0.499999);
for i=1:size(proj,3)
proj_CTF(:,:,i)=(proj(:,:,i)).*h(:, :, index(i));% generating data set with CTF (in Fourier space)
end

CTF=h;


