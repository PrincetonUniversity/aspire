classdef imagestackWriterComplex < imagestackWriter
% Out of core complex-valued image stack.
% Write images to a large image stack incrementally.
%
% For examples see test_imagestackWriter1, test_imagestackWriter2,
% test_imagestackWriter3.
%
% Yoel Shkolnisky, April 2017.


methods
    function obj = imagestackWriterComplex(filename,numslices,pixA,cachesize,mrcmode)
        % Constructor
        % filename  Name of MRC file to store the images on disk.
        % numslices Total number of slices that will be written to the
        %           stack. The current implementation of the MRC io
        %           functions requires to know this number in advance.
        % pixA      Pixel size of the images.
        % cachesize Number of images to store in memory.
        % mrcmode   Encoding used to the MRC file (0: uint8, 1: int16,
        %           2: single (default), 6: uint16).
        
        if ~exist('pixA','var')
            pixA=1;
        end
            
        if ~exist('cachesize','var') || isempty(cachesize)
            cachesize=1;
        end
        if ~exist('mrcmode','var')
            mrcmode=2;
        end
        
        obj=obj@imagestackWriter(filename,2*numslices,pixA,2*cachesize,mrcmode);
        obj.numslices=numslices;
    end
    
    % Append image to the stack.
    % The image is stored in memory. Once cache fills up, the images
    % are flushed to disk.
    function append(obj,images)
        % Write the given images into the MRC file (writing is delayed
        % until cache is full).
        
        imr=real(images);
        imc=imag(images);
        images=zeros(size(images,1),size(images,2),2*size(images,3));
        images(:,:,1:2:end)=imr;
        images(:,:,2:2:end)=imc;
        
        append@imagestackWriter(obj,images);
    end
    
end    
end
        