function [projs,rots]=cryo_project_random(vol,N)
% CRYO_PROJECT_RANDOM  Fast algorithm for generating randomly-oriented
%                      projections
%
% [projs,rots]=cryo_project_random(vol,N)
%   Generate randomly oriented projection of vol. The number of generated
%   projections is the nearest larger multiple of 3 of N (ceil(N/3)*3)). 
%   The function rotations vol in ceil(N/3) random orientations, for each
%   each orientation sums the resulting rotated volume in the X, Y, and Z
%   directions.
%
% Example:
%   load cleanrib
%   vol=real(volref); vol=GaussFilt(vol,0.5);
%   N=100;
%   [projs,rots]=cryo_project_random(vol,N);
%   projs2=cryo_project(vol,rots);
%   projs2=permute(projs2,[2,1,3]);
%   err=norm(projs(:)-projs2(:))/norm(projs(:))
%
%   Error should be about 1.0e-2.
%   Note that this function is not as accurate as cryo_project, but is
%   sufficient for most purposes.
%
% Yoel Shkolnisky, June 2019.

projs=zeros(size(vol,1),size(vol,2),N);
pv=permute(vol,[2 1 3]); % Permuted volume

rots_base = rand_rots(ceil(N/3));
roundedN=ceil(N/3)*3; % Round the number of projections to the
% closest multiple of 3, due to the way we generate projections.
rots=zeros(3,3,roundedN);

idx=1; % Current reference rotation we are generating
for k=1:ceil(roundedN/3)
    pv1= fastrotate3d(pv,rots_base(:,:,k));
    
    % Project in the Z direction
    projs(:,:,idx)=sum(pv1,3);
    rots(:,:,idx)=rots_base(:,:,k);
    idx=idx+1;
    % Compare the two images
    % imagesc(projs2(:,:,idx-1)); axis image; colormap(gray);    
    % imagesc(sum(fastrotate3d(pv,rots_ref(:,:,idx-1)),3)); axis image; colormap(gray);    


    % project in the X direction
    projs(:,:,idx)=squeeze(sum(pv1,1));
    rots(:,:,idx)=rotz(90)*rotx(90)*rots_base(:,:,k);
    idx=idx+1;
    % Compare the two images
    % imagesc(projs2(:,:,idx-1)); axis image; colormap(gray);
    % imagesc(sum(fastrotate3d(pv,rots_ref(:,:,idx-1)),3)); axis image; colormap(gray);    

    % Project in the Y direcition
    projs(:,:,idx)=squeeze(sum(pv1,2));
    rots(:,:,idx)=roty(90)*rots_base(:,:,k);
    idx=idx+1;
    % imagesc(projs2(:,:,idx-1)); axis image; colormap(gray);
    % imagesc(sum(fastrotate3d(pv,rots_ref(:,:,idx-1)),3)); axis image; colormap(gray);    
   
end
