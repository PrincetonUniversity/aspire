classdef imagestackReaderComplex < imagestackReader
% Out of core complex-calued image stack.
% Read complex-valued images from an file on demand.
% At any given time only cachesize images are stored in memory (default
% 100 complex arrays).
%
% Example:
%       stack=imagestackReaderComplex('FFTdata.mrc');
%       projs=zeros(65,65,100);
%       for k=1:100;
%            projs(:,:,j)=stack.getImage(k);
%       end
%
% Yoel Shkolnisky, April 2017.


methods
    function obj = imagestackReaderComplex(MRCname,cachesize,precision,verbose)
        % Construct and image stack object corresponding to a given
        % MRC filename.
        % cachesize     Number of images to store in memory
        % precision     'single' (default) or' double.
        % verbose       nonzero to print messages (default 0).
        
        if ~exist('cachesize','var')
            cachesize=100; % Default cache size.
        end
        
        if ~exist('precision','var')
            precision='single';
        end
        
        if ~exist(MRCname,'file')
            error('%s does not exist',MRCname);
        end
        
        if ~exist('verbose','var')
            verbose=0;
        end
                
        obj=obj@imagestackReader(MRCname,2*cachesize,precision,verbose);
        obj.dim(3)=obj.dim(3)/2;
    end
    
    function im=getImage(obj,idx)
        % Get image with index idx from the image stack.
        % If image is not in stack, flush entire stack and read images
        % from disk, starting with image idx. idx can be a vector.
        im=zeros(obj.dim(1),obj.dim(2),numel(idx),obj.precision);
        
        for k=1:numel(idx)
            imr=getImage@imagestackReader(obj,2*idx(k)-1);
            imc=getImage@imagestackReader(obj,2*idx(k));
            im(:,:,k)=imr+1i.*imc;
        end
    end
end
end