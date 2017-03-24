function rmin=cryo_masking_radius_3d_ks(V,verbose)
%
% This function does not work for reconstrudted density maps.
% Probably due to the properties of the noise in the reconstruction.
%
% CRYO_MASKING_RADIUS_KS   Estimate masking radius of a projection.
%
% rmin=cryo_masking_radius_ks(P)
%   Automatically estimate the masking radius of the volume V.
%   The function searches for the smallest centered sphere such that all
%   samples outside the sphere are noise pixels. This sphere is the
%   bounding circle of the molecule.
%   The function returns the radius of the bounding sphere.
%   The estimate is based on the Kolmogorov-Smirnoff test.
%
% rmin=cryo_masking_radius_ks(P,verbose)
%   Plot progress/debug figures. Default is verbose=0.
%
% Yoel Shkolnisky, October 2016.

if ~exist('verbose','var')
    verbose=0;
end

if ndims(V)~=3
    error('Input must be a 3D square array');
end
    
p=size(V,1);
if p~=size(V,2) || p~=size(V,3)
    error('Input must be a 3D square array');
end

% All pixels outside a circle of diameter p (pixels in the "corners" of the
% image) are considered noise pixels. Extract those pixels and use them as
% reference noise samples.
noise_radius = floor(p/2)-1; % Pixels outside this radius are considered 
                             % noise pixels
c=(p+1)/2; % Center of the image
[xmap,ymap,zmap]=ndgrid(1:p,1:p,1:p); % The i'th pixel in the volume is 
    % V(xmap(i),ymap(i),zmap(i)). xmap, ymap, and zmap are used to easily
    % find the pixels whose 3D coordinates satisfy some condition. See next
    % line.
noise_idx=find((xmap-c).^2+(ymap-c).^2+(zmap-c).^2>=noise_radius.^2); 
    % Find all "noise pixels", that is, all pixels outside the given radius.
X=V(noise_idx);

alpha=0.05; % All subsequent hypothesis tests will use this significance 
    % level. Smaller values of alpha rejects the null hypothesis less
    % often, so H==0 (see below) more often, and we get a tighter bounding
    % circle.
    
    % Make sure the image is not all noise
    Y=V(setdiff(1:p^3,noise_idx));
    H=kstest2(X,Y,'Alpha',alpha);
    
    if ~H
        log_message('Image is all noise.');
        rmin=-1;
        return;
    end
    
    rmin=p; % The smallest bouding radius detected.
    
    % Uncomment the following four lines if you want to search for a
    % non-centered bounding circle. This runs slower. Also remember to
    % return the detected center cx, cy, and cz.
    
    % cstep=ceil(p*0.05);   % Step used in center search.
    % for cntrx=floor(p/3):cstep:ceil(2*p/3)
    %    for cntry=floor(p/3):cstep:ceil(2*p/3)
    %       for cntrz=floor(p/3):cstep:ceil(2*p/3)
    for cntrx=c
        for cntry=c
            for cntrz=c
                distsq1=(xmap-cntrx).^2+(ymap-cntry).^2+(zmap-cntrz).^2;
                
                % Determine a circle large enough to enclose the central part of
                % the image.
                rmax=max([cntrx cntry cntrz p-cntrx p-cntry p-cntrz]);
                
                % Binary search of the masking raduis
                a=round(p*0.2); % Don't allow the circle to be too small.
                b=rmax;
                r=a; % Dummy value to be overriden next.
                
                while b-a>1
                    r=round((a+b)/2);
                    
                    % Take all samples in a shell of inner radius r and
                    % outer radius r+d, and check if they are distributed
                    % like the noise samples. We are looking to the
                    % smallest such ring (smallest r).
                    d=round(p*0.1); % Width of the ring is 10% the side of the image.
                    samples_idx=find(distsq1>=r.^2 & distsq1<=(r+d).^2);
                    
                    samples_idx=setdiff(samples_idx,noise_idx); % Remove any
                    % reference noise pixels from this set.
                    
                    if isempty(samples_idx)
                        % If no pixels left in the set, then probably the ring is
                        % too large, so set H=0 to make it smaller in the next
                        % iteration.
                        H=0;
                    else
                        Y=V(samples_idx);
                        [H,pval]=kstest2(X,Y,'Alpha',alpha);
                    end
                    
                    if H==0
                        % If the pixels in the ring are indistiguishable from
                        % noise, that is, they contain no parts of the object, then
                        % make the ring smaller and check again.
                        b=r;
                    else
                        % The samples in the ring can be distiguished from the
                        % reference noise samples. This may mean that they contain
                        % object pixels, that is, the ring should be enlarged to
                        % exclude object pixels.
                        a=r;
                    end
                    
                    if verbose
                        cryo_masking_radius_3d_ks_auxplot(V,H,r,[cntrx cntry cntrz]);
                        log_message('r=%d  H=%d  Pval=%d',r,H,pval);
                    end
                end
                
                % If we found a smaller circle, then update rmin.
                if r<rmin
                    rmin=r;
                    %cx=cntrx; % Uncomment for center search
                    %cy=cntry; % Uncomment for center search
                end
            end
        end
    end

cx=c; cy=c; cz=c;
if verbose
    cryo_masking_radius_3d_ks_auxplot(V,0,rmin,[cx cy cz])
end