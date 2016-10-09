function rmin=cryo_masking_radius_2d(P,verbose)
%
% CRYO_MASKING_RADIUS   Estimate masking radius of a projection.
%
% rmin=cryo_masking_radius_2d(P)
%   Automatically estimate the masking radius of the projection image P.
%   The function searches for the smallest centered circle such that all
%   samples outside the circle are noise pixels. This circle is the
%   bounding circle of the molecule.
%   The function returns the radius of the bounding circle.
%
% rmin=cryo_masking_radius_2d(P,verbose)
%   Plot progress/debug figures. Default is verbose=0.
%
% Yoel Shkolnisky, October 2016.

if ~exist('verbose','var')
    verbose=0;
end

if ~ismatrix(P)
    error('Input must be a 2D square image');
end
    
p=size(P,1);
if p~=size(P,2)
    error('Input must be a 2D square image');
end

% All pixels outside a circle of diameter p (pixels in the "corners" of the
% image) are considered noise pixels. Extract those pixels and use them as
% reference noise samples.
noise_radius = floor(p/2)-1; % Pixels outside this radius are considered 
                             % noise pixels
c=(p+1)/2; % Center of the image
[rmap,cmap]=ndgrid(1:p,1:p); % The i'th pixel in the image is 
    % P(rmap(i),cmap(i)). rmap and cmap are used to easily find the pixels
    % whose 2D coordinates satisfy some condition. See next line.
noise_idx=find((rmap-c).^2+(cmap-c).^2>=noise_radius.^2); % Find all "noise
    % pixels", that is, all pixels outside the given radius.
X=P(noise_idx);

% % Plot the pixels used as noise in white.
% maskim=zeros(p);
% maskim(noise_idx)=1;
% imagesc(maskim);
% colormap(gray);
% axis image;

alpha=0.05; % All subsequent hypothesis tests will use this significance 
    % level. Smaller values of alpha rejects the null hypothesis less
    % often, so H==0 (see below) more often, and we get a tighter bounding
    % circle.
    
% Make sure the image is not all noise
Y=P(setdiff(1:p*p,noise_idx));
H=kstest2(X,Y,'Alpha',alpha);

if ~H
    log_message('Image is all noise.');
    rmin=-1;
    return;
end

rmin=p; % The smallest bouding radius detected.

% Uncomment the following three lines if you want to search for a
% non-centered bounding circle. This runs a bit slower. Also remember to
% return the detected center cx and cy.

% cstep=ceil(p*0.05);   % Step used in center search.
% for cntrx=floor(p/3):cstep:ceil(2*p/3)
%    for cntry=floor(p/3):cstep:ceil(2*p/3)
for cntrx=c
    for cntry=c
        distsq1=(rmap-cntrx).^2+(cmap-cntry).^2;
        
        % Determine a circle large enough to enclose the central part of
        % the image.
        rmax=max([cntrx cntry p-cntrx p-cntry]);
        
        % Binary search of the masking raduis
        a=round(p*0.2); % Don't allow the circle to be too small.
        b=rmax;
        r=a; % Dummy value to be overriden next.
        
        while b-a>1            
            r=round((a+b)/2);

            % Take all samples in a ring of inner radius r and outer 
            % radius r+d, and check if they are distributed like the noise
            % samples. We are looking to the smallest such ring (smallest
            % r).
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
                Y=P(samples_idx);
                H=kstest2(X,Y,'Alpha',alpha);
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
                % Plot the currently tested bounding circle.
                imagesc(P);
                axis image;
                hold on;
                th = 0:pi/50:2*pi;
                xunit = r * cos(th) + cntrx;
                yunit = r * sin(th) + cntry;
                
                % Don't plot parts of the circle that fall outside the image.
                xunit(xunit<=0)=1; xunit(xunit>p)=p;
                yunit(yunit<=0)=1; yunit(yunit>p)=p;
                
                % If H==0, the our hypothesis was not rejected, that is, we
                % think that the ring contains only noise pixels, and so our
                % object is enclosed in it. In such a case the circle is green.
                % Otherwise, be think the circle intersects the object, and so
                % the circle is red.
                if H==0
                    clr='g';
                else
                    clr='r';
                end
                
                plot(xunit, yunit,'LineWidth',2,'Color',clr);
                title(sprintf('radius=%d',r));
                hold off;
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

cx=c;
cy=x;
if verbose
    % Plot the final detected circle
    imagesc(P);
    axis image;
    hold on;
    th = 0:pi/50:2*pi;
    xunit = rmin * cos(th) + cx;
    yunit = rmin * sin(th) + cy;
    
    % Don't plot parts of the circle that fall outside the image.
    xunit(xunit<=0)=1; xunit(xunit>p)=p;
    yunit(yunit<=0)=1; yunit(yunit>p)=p;
    
    plot(xunit, yunit,'LineWidth',2,'Color','g');
    hold off;
end