function [bestR,bestdx,bestCorr,bestResA,vol2aligned]=bf3Dmatchaux(vol1,vol2,rotations,pixA,verbose)
%
% BF3DMATCHAUX Auxiliary function for matching 2 density maps.
%
% This function should be called directly.
%
% [bestR,bestdx,bestCorr]=bf3Dmatchaux(vol1,vol2,rotations)
%   Use the given "rotations" to try to best align vol1 and vol2.
%   The function rotates vol2 by all rotations in the list "rotations", and
%   for each rotation finds the best translation. It returns the best
%   rotation fro the list, the corresponding best translation parameters,
%   and the correlation of the two volumes after the volumes haved been
%   aligned using these rotation and translation.
%
%   Input parameters:
%       vol1        nxnxn density map.
%       vol2        nxnxn density map.
%       rotation    List of rotations, where rotations(:,:,k) is the k'th
%                   rotation.
%
%   Output parameters:
%       bestR       Optimal rotation from the list "rotations".
%       bestdx      Optimal estimated translation.
%       bestCorr    Correlation between vol1 and vol2 after aligning them
%                   using bestR and bestdx.
%       vol2aligned vol2 after aligning it to vol1.
%
% Yoel Shkolnisky, January 2015.

if ~exist('verbose','var')
    verbose=0;
end

if verbose
    printProgressBarHeader;
end

fakepixA=0;
if pixA<=0
    pixA=1;      % Use dummy pixel size.
    fakepixA=1 ; % Indicate that we are using fake pixel size
end

%ResVec=zeros(size(rotations,3),1); % Resolution achieved by each rotation
%Nvec=zeros(size(rotations,3),1); % Number of elements in each FSC curve.


parfor jj=1:size(rotations,3)
%for jj=1:size(rotations,3)
    progressTic(jj,size(rotations,3));
    rot=rotations(:,:,jj);    
    vol2R=fastrotate3d(vol2,rot);    
    estdx=register_translations_3d(vol1,vol2R);
    if ~isscalar(estdx)
        vol2RS=reshift_vol(vol2R,estdx);
%         fsc=FSCorr(vol1,vol2RS);
%         res=fscres(fsc,0.134);
%         ResVec(jj)=res;
%         Nvec(jj)=numel(fsc); % To compute resolution in angstrom.
        Cvec(jj)=corr(vol1(:),vol2RS(:));
    end

end

[bestCorr,ii]=max(Cvec);
% [bestRes,ii]=max(ResVec);
% bestResA=2*pixA*Nvec(ii)/bestRes; % Resolution in Angstrom.
% if fakepixA
%     bestResA=-bestResA; % Minus is an indication of fake pixel size
% end
bestR=rotations(:,:,ii);
vol2aligned=fastrotate3d(vol2,bestR);
bestdx=register_translations_3d(vol1,vol2aligned);
vol2aligned=reshift_vol(vol2aligned,bestdx);    
%bestCorr=corr(vol1(:),vol2aligned(:));
fsc=FSCorr(vol1,vol2aligned);
res=fscres(fsc,0.134);
bestResA=2*pixA*numel(fsc)/res; % Resolution in Angstrom.
if fakepixA
    bestResA=-bestResA; % Minus is an indication of fake pixel size
end
