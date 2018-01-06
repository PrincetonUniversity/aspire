function err=cryo_compare_stacks(mrcname1,mrcname2,verbose)
%
% CRYO_COMPARE_STACKS   Compare between two MRC stacks.
%
% err=cryo_compare_stacks(mrcname1,mrcname2)
%   Compare two MRC stacks stored in the files mrcname1 and mrcname2 image
%   by image. Returns the relative error between the stacks
%
%   verbose=0   silent
%   verbose=1   show progress bar
%   verbose=2   print progress every 1000 images
%   verbose=3   print message for each processed image
%
% Yoel Shkolnisky, October 2017

if ~exist('verbose','var')
    verbose=0;
end

mrc1reader=imagestackReader(mrcname1);
mrc2reader=imagestackReader(mrcname2);

% Check that the dimensions of the stack are compatible
if mrc1reader.dim(1)~=mrc2reader.dim(1)
    error('x dimension in both stacks is not compatible: %s has %d pixels, %s is %d pixels',...
        mrcname1,mrc1reader.dim(1),mrcname2,mrc2reader.dim(1));
end

if mrc1reader.dim(2)~=mrc2reader.dim(2)
    error('y dimension in both stacks is not compatible: %s has %d pixels, %s is %d pixels',...
        mrcname1,mrc1reader.dim(2),mrcname2,mrc2reader.dim(2));
end

if mrc1reader.dim(3)~=mrc2reader.dim(3)
    error('z dimension in both stacks is not compatible: %s has %d pixels, %s is %d pixels',...
        mrcname1,mrc1reader.dim(3),mrcname2,mrc2reader.dim(3));
end

errtot=0;

if verbose==1
    printProgressBarHeader;
end

for k=1:mrc1reader.dim(3)
    if verbose==1
        progressTicFor(k,mrc1reader.dim(3));
    end

    im1=mrc1reader.getImage(k);
    im2=mrc2reader.getImage(k);
    
    ek=norm(im1(:)-im2(:));
    
    if verbose>=3
        log_message('Difference between projections %d/%d: %5.3e',ek);
    end

    errtot=errtot+ek;

    if verbose==2
        if mod(k,1000)==0
            log_message('Finsihed processing %d/%d projections. Relative error so far %5.3e'...
                ,k,mrc1reader.dim(3),errtot/k);
        end
    end    
end

err=errtot/mrc1reader.dim(3);
if verbose>=3
    log_message('Relative error between stacks: %5.3e',err);
end

