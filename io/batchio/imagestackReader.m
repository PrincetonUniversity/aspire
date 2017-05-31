classdef imagestackReader < handle
    % Out of core image stack.
    % Read images from an MRC file on demand.
    % At any given time only cachesize images are stored in memory (default
    % 100).
    %
    % Example:
    %       stack=imagestackReader('stack.mrc');
    %       projs=zeros(65,65,100);
    %       for k=1:100;
    %            projs(:,:,j)=stack.getImage(k);
    %       end
    %
    % Yoel Shkolnisky, May 2016.
    
    properties
        filename            % Name of the filename containing the image.
        images              % Cache of images loaded into memory.
        images_idx          % Indices of the loaded images.
        cachesize           % Number of images to store in memory.
        precision           % Precision of the stored data. Also detetmines 
                            % the type of the returned data.
        dim                 % Dimensions of the image stack. First two
        % dimensions are image size. Third dimension is number of images.
        % These are the logical dimensions. For example, for N complex
        % images, the third dimension of dim would be N, but the third
        % dimension of physical_dim below would be 2N.
        % the 
        verbose             % Print mesages.
    end
    
    properties (Access=private)
        physical_dim    % Physical dimensions of the stored data.
    end
    
    methods
        function obj = imagestackReader(MRCname,cachesize,precision,verbose)
            % Construct an image stack object corresponding to a given
            % MRC filename.
            % cachesize     Number of images to store in memory
            % precision     'single' (default) or' double.
            % verbose       nonzero to print messages (default 0).
            
            if ~exist('cachesize','var')
                cachesize=100; % Default cache size.
            end
            obj.cachesize=cachesize;
            
            if ~exist('precision','var')
                precision='single';
            end
            
            if ~exist(MRCname,'file')
                error('%s does not exist',MRCname);
            end
            
            if ~exist('verbose','var')
                verbose=0;
            end
            
            [im,info]=ReadMRC(MRCname,1,1); % Read the first image;
            obj.filename=MRCname;
            obj.precision=precision;
            obj.images=zeros(size(im,1),size(im,2),obj.cachesize,obj.precision);
            obj.images_idx=zeros(obj.cachesize,1);
            obj.dim=[size(im,1);size(im,2);info.nz;];
            obj.physical_dim=obj.dim;
            obj.verbose=verbose;
            
            if obj.verbose
                log_message('Creating imagestack. Nimages=%d, Ncache=%d, precision=%s, MRC=%s',...
                    obj.dim(3),obj.cachesize,precision,obj.filename);
            end
            
        end
        
        function im=getImage(obj,idx)
            % Get image with index idx from the image stack.
            % If image is not in stack, flush entire stack and read images
            % from disk, starting with image idx. idx can be a vector.
            im=zeros(obj.physical_dim(1),obj.physical_dim(2),numel(idx),obj.precision);
            
            for k=1:numel(idx)                
                if idx(k)>obj.physical_dim(3)
                    error('Image %d does not exist in stack %s',idx(k), obj.filename);
                end
                
                ii=find(obj.images_idx==idx(k));
                if isempty(ii)
                    % Read images
                    if obj.verbose
                        log_message('Cache miss for image %d. Loading images.',idx(k));
                    end
                    
                    buf=ReadMRC(obj.filename,idx(k),obj.cachesize);
                    nread=size(buf,3);
                    obj.images(:,:,1:nread)=buf;
                    obj.images_idx=idx:idx+nread-1;
                    ii=1;
                end
                im(:,:,k)=obj.images(:,:,ii);
            end
        end
    end
end
