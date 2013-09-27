function [m olmap]=RemoveOutliers(m,nsd,maxfract,showWarnings)
% use: [out olmap]=RemoveOutliers(m,nsd,maxfract,showWarnings);
% or:        out = RemoveOutliers(m,olmap);
% Remove outlier pixels from images by forming local medians of
% neighboring pixels. nsd is the number of standard deviations at which
% outliers are defined; default is 5.  maxfract is the fraction of
% points allowed to be outliers (default is 10^-4).  The more stringent of
% these two criteria is used, so to just select a fraction you can set
% nsd=0. The last argument showWarnings (default 0) turns on a warning
% when nearby valid pixels are not found.
% The returned olmap is 1 wherever a pixel was corrected.
% To correct a second image based on the pattern, use the second form,
% passing the previous olmap.
% 
% If m is an integer array, the values are converted to single.
% fs 2 sep 07, revised 30 Jun 09, 7 Sep 10.
% revised to take medians and use Percentile() 26 May 12.

% Search Parameters
minStep=2;  % start with a 5x5 region
maxStep=8;  % maximum size of neighborhood to search is (2*maxstep+1) square
stepIncrement=2;
minPoints=3;   % Get the median of at least this many points


msize=size(m);
nx=msize(1);
ny=msize(2);
if isinteger(m)  % convert to floating-point
    m=single(m);
end;
if nargin<4
    showWarnings=0;
end;
if nargin<3
    maxfract=1e-4;
end;
if nargin<2
    nsd=5;  % number of standard deviations for allowed points.
end;

me=mean(m(:));
if numel(nsd)<2  % a scalar; we define the outliers.
    sd=std(m(:));
    [lowerLimit upperLimit]=Percentile(m,maxfract);
    upperLimit=max(me+nsd*sd,upperLimit);
    lowerLimit=min(me-nsd*sd,lowerLimit);
    olmap=(m>upperLimit)|(m<lowerLimit);  % find all the values outside the limits.
else
    olmap=nsd;  % alternative use of the variable is the outlier map.
end;

% we now have the outlier map find the nonzero elements.
ol=find(olmap);

% Pick up the subscripts of the outlying points
[is js]=ind2sub(msize,ol);

errno=0;    % number of failures

for k=1:numel(ol);
    sn=0;
    istep=minStep;
    i=is(k); j=js(k);  % Coordinate of the point in the original image
    while (sn<minPoints) && (istep<=maxStep) % search a neighborhood
        % make a copy of the points surrounding the outlier
        imin=max(1,i-istep);
        imax=min(nx,i+istep);
        jmin=max(1,j-istep);
        jmax=min(ny,j+istep);
        mlocal=m(imin:imax,jmin:jmax);
        goodlocal=(olmap(imin:imax,jmin:jmax)==0);
        sn=sum(goodlocal(:));
        if sn>=minPoints
            s=median(mlocal(goodlocal(:)>0));
        end;
        istep=istep+stepIncrement;              % neighborhood size
    end;
    % Replace the outlier with the median of in-bounds surrounding points.
    if sn<minPoints  % no neighborhood found, use the global mean.
        errno=errno+1;
        m(i,j)=me;
    else
        m(i,j)=s;
    end;
    
end;
if showWarnings && errno>0
    warning('Neighborhood couldn''t be found for %d points',errno);
end;

