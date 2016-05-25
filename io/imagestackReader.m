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
        stacksize           % Total number of images in the stack.
        verbose             % Print mesages.
    end
        
    methods
        function obj = imagestackReader(MRCname,cachesize,precision,verbose)
            % Construct and image stack object corresponding to a given
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
            obj.images=zeros(size(im,1),size(im,2),obj.cachesize,precision);
            obj.images_idx=zeros(obj.cachesize,1);   
            obj.stacksize=info.nz;
            obj.verbose=verbose;
            
            if obj.verbose
                log_message('Creating imagestack. Nimages=%d, Ncache=%d, precision=%s, MRC=%s',...
                    obj.stacksize,obj.cachesize,precision,obj.filename);
            end
                
        end
    
        function im=getImage(obj,idx)
            % Get image with index idx from the image stack.
            % If image is not in stack, flush entire stack and read images
            % from disk, starting with image idx.
            
            if idx>obj.stacksize
                error('Image %d does not exist in stack.',idx);
            end
            
            ii=find(obj.images_idx==idx);
            if isempty(ii)
                % Read images
                if obj.verbose
                    log_message('Cache miss. Loading images.');
                end
                
                buf=ReadMRC(obj.filename,idx,obj.cachesize);
                nread=size(buf,3);
                obj.images(:,:,1:nread)=buf;
                obj.images_idx=idx:idx+nread-1;
                ii=1;
            end
            im=obj.images(:,:,ii);            
        end
    end
end