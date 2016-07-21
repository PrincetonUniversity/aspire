classdef gpumatrix < handle
% GPU varaiable.
%
% A pointer to a matrix stored on the GPU.
%
% Yoel Shkolnisky, July 2016.

properties
        gptr                % Pointer into the GPU memory.
        dim                 % Dimensions of the matrix.
        mode                % 0 for single, 1 for single-complex.
    end
        
    methods
        function obj = gpumatrix(varargin)
            % Create a gpumatrix object.
            %   gpumatrix(A) Create a gpumatrix object from a matrix A.
            %       Store the matrix A on the GPU and save the handle to
            %       the stored GPU object. The matrix A must be single
            %       precision real or complex. 
            %   gpumatrix(gptr,M,N,mode) Create a gpumatrix object from a
            %       deivce pointer 
            
            if nargin==1
                A=varargin{1};
                if ~ismatrix(A)
                    error('Only matrices are supported');
                end
                
                if ~isa(A,'single')
                    error('Only single precision is supported');
                end
                
                obj.mode=0;
                if ~isreal(A)
                    obj.mode=1;
                end
                
                obj.dim=size(A);
                
                if obj.mode==0
                    obj.gptr=gpuauxSalloc(A);
                else
                    obj.gptr=gpuauxCalloc(A);
                end
            elseif nargin==4
                obj.gptr=varargin{1};
                obj.dim(1)=varargin{2};
                obj.dim(2)=varargin{3};
                obj.mode=varargin{4};
            else
                error('Incorrect number of arguments.');
            end
        end
        
        function delete(obj)
            % Destructor. Free resources for the GPU object.
            free(obj);
        end

        function free(obj)
            % Free the stored matrix from the GPU and invalidate the
            % object to prevent further access.
            if obj.gptr~=0
                gpuauxfree(obj.gptr);
                obj.gptr=0;
            end
        end
        
        function A=gather(obj) 
            % Read the matrix from the GPU.
            if obj.gptr==0
                error('No GPU data exists');
            end
            
            if obj.mode==0
                A=gpuauxSgather(obj.gptr,obj.dim(1),obj.dim(2));
            else
                A=gpuauxCgather(obj.gptr,obj.dim(1),obj.dim(2));
            end
        end
    end
    
end