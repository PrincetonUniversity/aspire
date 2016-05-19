function [estR,est_shifts,rawimageindices,stats]=cryo_assign_orientations_to_raw_projections(...
    alignment_data,averaging_data,ref_R,ref_shifts)
%
% CRYO_ASSIGN_ORIENTATIONS_TO_RAW_PROJCETIONS Assign orientations to the
%   individual images comprising the class averages.
% [estR,est_shifts,rawimageindices,stats]=cryo_assign_orientations_to_raw_projections(...
%       alignment_data,averaging_data)
%   Given alignment data returned from the class averaing procedure
%   together with the averages and their oreintation parameters, compute
%   for each indivial image in each of the class averages its orientation
%   and shift.
%
% [estR,est_shifts,rawimageindices,stats]=cryo_assign_orientations_to_raw_projections(...
%       alignment_data,averaging_data,ref_R,ref_shifts)
%   If reference orientations and reference shifts are given, they are used
%   to compute the errors in the estimated rotations and shifts.
%
% Input parameters:
%   alignment_data  Structure with the fields classcoreidx, class_VDM,
%                   VDM_angles, shifts, class_VDM_refl. There information
%                   is generated during the VDM classification process.
%   averaging_data  Structure with the fields rotations, est_shifts, nnavg.
%                   That is, the rotation and estimated shift for each
%                   class averages. Each average was generated from nnavg
%                   raw images.
%   ref_R           Reference rotation of each raw image (for debugging).
%   ref_shift       Reference shift of rach raw image (for debugging).
%
% Output parameters:
%   estR                Estimated rotations of each raw image.
%   est_shifts          Estimated shift of each raw image.
%   rawimageindices     Indices of raw images that participate in the the
%                       given averages. Rotations and shifts are estimated
%                       only for those images.
%   stats               Statistics on the assignment process. Structure
%                       with the fields rawimage_num_count, Rs_spectrum,
%                       shifts_deviations.
% 
% Structure stats:
%   rawimage_num_count  Number of times each raw image was used. Array of
%                       length numel(rawimageindices) where
%                       rawimage_num_count(i) is the number of times image
%                       rawimageindices(i) was used for averaging.
%   Rs_spectrum         Singular values of the mean rotation for raw images
%                       that participate in several class averages and so
%                       have several estimated rotations.
%   shifts_deviations   Deviation of the estimated shifts from the median
%                       shift for raw images that participate in several
%                       class averages.
%
% Yoel Shkolnisky, May 2016.
nmeans=size(averaging_data.rotations,3);
rotations=averaging_data.rotations(:,:,1:nmeans);
nnavg=averaging_data.nnavg;

ref_R_given=0;
if exist('ref_R','var')
    ref_R_given=1;
    % If reference rotations are given, regsiter the estimted rotations to the
    % reference ones. Note that we have to reorder the reference rotations in
    % the same order as the averages, namely, according to
    % alignment_data.classcoreidx. This is why ref_R is reordered as
    % ref_R(:,:,alignment_data.classcoreidx(1:nmeans)).
    [rotations,mse,~,~,Jflag]=register_rotations(rotations,...
        ref_R(:,:,alignment_data.classcoreidx(1:nmeans)));
    % mse should be very small.
    % flag is 2 is there is J-conjugation between the two sets of rotations and
    % 1 otherwise.
    log_message('Reference rotations are given. Registering reference rotations to given ones.');
    log_message('Registration: mse=%d, J-conjugation flag=%d',mse,Jflag)
end

ref_shifts_given=0;
if exist('ref_shifts','var')
    ref_shifts_given=1;
end

% Each row in alignment_data.class_VDM contains the nearest neighbors of
% the image corresponding to the row. alignment_data.classcoreidx(1:nmeans)
% are the indices of the images that used to generate the first nmeans
% averages.
% alignment_data.class_VDM(alignment_data.classcoreidx(1:nmeans),:) are the
% nearest neighbors of these images.
NNidx=alignment_data.class_VDM(alignment_data.classcoreidx(1:nmeans),:);
NNidx=NNidx(:,1:nnavg);

% Find the indices of the raw images used to generate the class averages.
rawimageindices=unique(NNidx); % These are the indices of the raw 
            % projections that should be used for reconstruction together
            % with the estimated rotations and shifts.
Nrawimageidx=numel(rawimageindices);
log_message('Class averages were generated from %d raw projections.',Nrawimageidx);

estR=zeros(3,3,Nrawimageidx);
est_shifts=zeros(Nrawimageidx,2);
%ref_shifts_subset=zeros(Nrawimageidx,2);
rawimage_num_count=zeros(Nrawimageidx,2); % Number of times each raw image 
                                          % was used.
Rs_spectrum=zeros(Nrawimageidx,3); % Singular values of the mean rotation 
    % for raw images that participate in several class averages and so have
    % several estimated rotations.
shifts_deviations=zeros(Nrawimageidx,2); % For raw images that participate 
    % in several class averages, we compute the deviation of each estimated
    % shift from the median shift.
% For each image that takes part in at least one of the averages, estimate
% its rotation based on its in-plane rotation with respect to the
% average(s) together with the estimated 3D rotation of the average(s).
% This results in several estimated of the 3D rotation of the image.

for k=1:Nrawimageidx; % For each of the raw images used...    
    currrawimage=rawimageindices(k); % Find the index of the next image that participates in at least one average.
    idx=find(NNidx==currrawimage);   % Find all averages in which it takes part.
    
    log_message ('Raw image %d appears %d times',currrawimage,numel(idx));
    rawimage_num_count(k,1)=currrawimage; % Index of the image
    rawimage_num_count(k,2)=numel(idx); % How many times it appeared.
    
    [I,J]=ind2sub(size(NNidx),idx); % The image whose index in the original 
        % stack of images is currrawimage participates in the averages
        % whose indices are given by the list I. The indices of the cores
        % of these averages (the core is the raw image from which the class
        % average was generated by averging it with its nearest neighbors)
        % in the original stack of images is
        % alignment_data.classcoreidx(I(:)). Note that
        % NNidx(I(l),J(l))=currrawimage. 
    
    % Compute the rotation and shift of each occurence of currrawimage in
    % NNidx.
    Rs=zeros(3,3,numel(idx)); % Estimates of the rotations of currrawimage.
    ds=zeros(numel(idx),2);   % Estimates of the shifts of currrawimage. 
    ref_ds=zeros(numel(idx),2); % Reference shifts of all occurences of 
        %currimageidx. Note that if dsj contains more than on row, then all
        %rows are identical, since all correspond to the same raw image
        % currrawimage.

    normRerrsq=0; % For debug - the difference between the estimated 
                 % rotations and the true one.
                 
    for l=1:numel(idx); % For each occurence of currrawimage...
        R1=rotations(:,:,I(l)); % Take the estimated 3D rotation of the 
            % average to which the current appearance of currrawimage was
            % aligned. 
        theta=alignment_data.VDM_angles(alignment_data.classcoreidx(I(l)),J(l));
            % Take the planar rotation angle between currrawimage and the
            % average.
        
        % For debug, we can use the following line that uses the reference
        % shifts.
        %ds1=ref_shifts(alignment_data.classcoreidx(I(l)),:); 
        %ds1=flipud(ds1(:));
        
        ds1=averaging_data.est_shifts(I(l),:); ds1=flipud(ds1(:));
            % Estimated shift for the class average to which the current
            % appearance of currimageidx was aligned (see R1 above).
        
        ii=alignment_data.classcoreidx(I(l)); jj=J(l);
        ds1l=[imag(alignment_data.shifts(ii,jj+1));...
            real(alignment_data.shifts(ii,jj+1)); 0];
            % Shift required to align the current appearance of
            % currrawimage to the its class avreage. Note that the y
            % coordinate given by the IMAG part goes before the x
            % coordinate given by the REAL part, in accordance with the
            % flipud in ds1.

                                     
        % Combine R1 and the in-plane rotation theta in an estimate of the
        % 3D rotation of the current appearance of currrawimage. We also
        % need to take care of the case where there is reflection between
        % currrawimage and the average.
        if alignment_data.class_VDM_refl(alignment_data.classcoreidx(I(l)),J(l))==1
            Rs(:,:,l)=R1*Rz(theta*pi/180);  % No reflection.
            tmp=Rz(-theta*pi/180)*([ds1;0]-ds1l);
            assert(tmp(3)==0)
            ds(l,:)=tmp(1:2).';
        else
            Rs(:,:,l)=R1*Rz(theta*pi/180)*diag([1,-1,-1]); % With reflection.
            tmp=diag([1,-1,-1])*Rz(-theta*pi/180)*([ds1;0]-ds1l);
            assert(tmp(3)==0)
            ds(l,:)=tmp(1:2).';
        end

        if ref_R_given
            % Compare to reference rotations
            % Compute the error of the estimated rotation relative to the true
            % one.
            Rii=ref_R(:,:,alignment_data.class_VDM(alignment_data.classcoreidx(I(l)),J(l)));
            normRerrsq=normRerrsq+norm(Rii.'*Rs(:,:,l)-eye(3)).^2;
        end
        
        if ref_shifts_given
            ii=alignment_data.classcoreidx(I(l)); jj=J(l);        
            ref_ds(l,:)=fliplr(ref_shifts(alignment_data.class_VDM(ii,jj),:));
        end
    end
    
    if ref_R_given
        e=sqrt(normRerrsq)/numel(idx);
        log_message('error in rotations=%d',e);
        %assert(e<1)
    end
    
    % Combine all estimated rotations of the different appearances of
    % currrawimage into a single estimate. The estimate R satisfies 
    %       min \sum_{l} ||R-Rs(l)||_{F}^{2}
    [U,S,V]=svd(mean(Rs,3));
    R=U*V.';
    estR(:,:,k)=R;
    log_message('Singular values of the mean of Rs = [%4.2e, %4.2e, %4.2e]',...
        S(1,1),S(2,2),S(3,3)); % If all estimates in Rs are very close to 
        % each other then all singular values should be very close to 1.
    Rs_spectrum(k,:)=diag(S);
% We cannot directly compare est_shifts and ref_est_eshifts as the two sets
% of shifts are not registered. To compare them we need the shift
% equations, as it done in cryo_estimate_shifts.
%     if norm(dsj,'fro')>0
%         e=norm(ds-dsj,'fro')/norm(dsj,'fro');
%     else
%         e=norm(ds-dsj,'fro');
%     end
%     log_message('error in shifts=%d',e);
%     assert(e<1)

if size(ds,1)>1
    est_shifts(k,:)=fliplr(median(ds));
    serr=sqrt(sum((abs(bsxfun(@minus,fliplr(ds),est_shifts(k,:)))).^2));
    log_message('norm of shift deviations = [%4.2e, %4.2e]',serr(1),serr(2));
    shifts_deviations(k,:)=serr;
else
    est_shifts(k,:)=fliplr(ds);
end

% ii=alignment_data.classcoreidx(I(1));
% ref_shifts_subset(k,:)=shifts(alignment_data.class_VDM(ii,1));

end

stats.rawimage_num_count=rawimage_num_count;
stats.Rs_spectrum=Rs_spectrum;
stats.shifts_deviations=shifts_deviations;