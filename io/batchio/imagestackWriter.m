classdef imagestackWriter < handle
% Out of core image stack.
% Write images to a large image stack incrementally.
%
% For examples see test_imagestackWriter1, test_imagestackWriter2,
% test_imagestackWriter3.
%
% Yoel Shkolnisky, May 2016.

    properties
        filename            % Name of the filename containing the image.
        fid                 % Handle to open MRC file.
        imagebuf            % Buffer of images not yet written.
        cachesize           % Number of images to store in memory.
        ncache              % Number of images in cache.
        dim                 % Dimension of each image.
        numslices           % Total number of images in stack.
        nwritten            % Number of images written so far.
        headerwritten       % Was a head already written to the MRC file.
        pixA                % Pixel size in Angstrom
        mrcmode             % MRC data type. Default is 2 (single precision);
    end
    
    properties (Access=private)
        physical_numslices  % Number of slices actually written to the 
            % file. Different from numslices if complex-valued data is
            % written.
    end
    
    methods
        function obj = imagestackWriter(filename,numslices,pixA,cachesize,mrcmode)
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
                        
            
            obj.filename=filename;
            obj.cachesize=cachesize;
            obj.ncache=0;
            obj.numslices=numslices;
            obj.physical_numslices=numslices;
            obj.nwritten=0;
            obj.headerwritten=0;
            obj.pixA=pixA;
            obj.mrcmode=mrcmode;
        end
        
        function delete(obj)
            close(obj);
        end
        
        % Append image to the stack.
        % The image is stored in memory. Once cache fills up, the images
        % are flushed to disk.
        function append(obj,images)
            % Write the given images into the MRC file (writing is delayed
            % until cache is full).
            
            if ~isreal(images)
                error('Images to write must be real-valued. For complex-valued images use imagestackWriterComplex');
            end
            
            % Verify dimensions of images.
            sz=size(images);
            nimages=1;
            if ndims(sz)>3 || ndims(sz)<2
                error('images should be an image or an image stack');
            end
            if ndims(images)==3
                nimages=sz(3);
                obj.dim=sz(1:2);
            else
                sz(3)=1;
            end
            
            if ~obj.headerwritten
                % Write header, update dim and header flag, allcoate imagebuf,
                org=-floor(sz/2);  % Default origin
                % Send the first image as it is used to extract image params.
                obj.fid=WriteMRCHeader(images(:,:,1),obj.pixA,...
                    obj.filename,obj.physical_numslices,org,obj.mrcmode);
                
                obj.headerwritten=1;
                obj.dim=sz;
                obj.imagebuf=zeros(obj.dim);
                obj.ncache=0;
                obj.nwritten=0;
            end
            % Check dimensions of images and add images to cache.
            % If cache is above chache size then flush.
            if sz(1)~=obj.dim(1) || sz(2)~=obj.dim(2)
                error('Size of images different from previously written images.');
            end
            
            for k=1:nimages
                if obj.ncache==obj.cachesize
                    flush(obj);
                end
                obj.ncache=obj.ncache+1;
                obj.imagebuf(:,:,obj.ncache)=images(:,:,k);
                obj.nwritten=obj.nwritten+1;
                if obj.nwritten>obj.physical_numslices
                    error('Trying to add more images than set for this stack.');
                end
                
            end
        end
                
        function flush(obj)
            % Flush unwritten images from cache to disk
            if obj.ncache>0
                % Write cache to disk
                fwrite(obj.fid,obj.imagebuf(:,:,1:obj.ncache),imagestackWriter.mrcmodestr(obj.mrcmode));
                obj.ncache=0; % Reset cache
            end
        end        
        
        
        function close(obj)
            % Close the stack.
            flush(obj);
            
            try
                if obj.fid>0
                    fclose(obj.fid);
                end
            catch
                % Using try...catch in case the file handle was
                % invalidated, due to, for example, debug.
            end
        end
        
        
    end
    
    methods(Static)
        function string=mrcmodestr(mode)
            % Convert MRC mode indentifier to string.
            switch mode
                case 0
                    string='uint8';
                case 1
                    string='int16';
                case 2
                    string='float32';
                case 6
                    string='uint16';
                otherwise
                    error(['Invalid data mode: ' num2str(mode)]);
            end;
        end               
    end
end
        