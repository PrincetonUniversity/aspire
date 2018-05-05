function [shifted_projections,ref_shifts] = cryo_addshifts(projections,shifts,max_shift,shift_step)
%
% Add shifts to projection images.
%
% Input parameters:
%   projections     3D stack of projections. The slice projections(:,:,k)
%                   is the k'th projections.
%   shifts          (Optional) A two column table with the 2D shift of each
%                   projection. Number of rows must be equal to the number
%                   of proejctions. If this parameter is not provided, pass
%                   an empty array. If provided, the following to
%                   parameters can be omitted.
%   max_shift       (Optional) Maximal random shift (in pixels) introduced
%                   to each projection. max_shift must be an integer and
%                   the resulting random shiftsare integers between
%                   -max_shift to +max_shift. Default max_shift is 0 (no
%                   shift).  
%   shift_step      (Optional) Resolution used to generate shifts.
%                   shift_step=1 allows for all integer shifts from
%                   -max_shift to max_shift (in both x and y directions).
%                   shift_step=2 generates shifts between -max_shift and
%                   max_shift with steps of 2 pixels, that is, shifts of 0
%                   pixels, 2 pixels, 4 pixels, and so on. Default
%                   shift_step is 0. 
%
% Output parameters:
%   shifted_projections  Shifted stack of images of the same size as the
%                       projections input array. 
%   ref_shifts          A two column table with the 2D shift introduced to
%                       each projections.
%
% Yoel Shkolnisky, September 2013.

K=size(projections,3);

if ~exist('max_shift','var')
    max_shift=0;
end

if ~exist('shift_step','var')
    shift_step=0;
end

if ~isempty(shifts)
    if size(shifts,1)~=K    % There must be one shift for each projections
        error('malformed "shifts". Must be of the same length as "projections"');
    end
    ref_shifts=shifts;
else
    ref_shifts=zeros(K,2);
    if shift_step>0        
        ref_shifts=round((rand(K,2)-1/2)*2*max_shift/shift_step)*shift_step;
    end

end

% Determine dimensions of each projection.
nx = size(projections,1);
ny = size(projections,2);

if mod(nx,2)==1
    rangex=-(nx-1)/2:(nx-1)/2;    
else
    rangex=-nx/2:nx/2-1; % Note that in the case of an even image, the 
                         % center is not at the middle. This can be easily
                         % fixed by switching to the appropriate FFT
                         % routines.
end

if mod(ny,2)==1
    rangey=-(ny-1)/2:(ny-1)/2;    
else
    rangey=-ny/2:ny/2-1;
end

[omega_x,omega_y]=ndgrid(rangex,rangey);
omega_x=-2*pi.*omega_x/nx; omega_y=-2*pi.*omega_y/ny;

shifted_projections=zeros(size(projections));

for k=1:K
    p=projections(:,:,k);
    pf=fftshift(fft2(ifftshift(p)));
    phase_x=omega_x.*ref_shifts(k,1);
    phase_y=omega_y.*ref_shifts(k,2);
    pf=pf.*exp(sqrt(-1)*(phase_x+phase_y));
    p2=fftshift(ifft2(ifftshift(pf)));    
    shifted_projections(:,:,k)=real(p2);
end

