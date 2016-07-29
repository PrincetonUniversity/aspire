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
        
        function gc = mtimes(ga, gb)
            if isnumeric(gb) 
                c = single(gb);
                if ga.mode == 1
                    gcptr=gpuauxCconstMul(ga.gptr,ga.dim(1),ga.dim(2), c);
                    gc=gpumatrix(gcptr,gb.dim(1),gb.dim(2),1);
                else
                    gcptr=gpuauxSconstMul(ga.gptr,ga.dim(1),ga.dim(2), c);
                    gc=gpumatrix(gcptr,gb.dim(1),gb.dim(2),0);
                end
                return
            elseif isnumeric(ga)
                c = single(ga);
                if gb.mode == 1
                    gcptr=gpuauxCconstMul(gb.gptr,gb.dim(1),gb.dim(2), c);
                    gc=gpumatrix(gcptr,gb.dim(1),gb.dim(2),1);
                else
                    gcptr=gpuauxSconstMul(gb.gptr,gb.dim(1),gb.dim(2), c);
                    gc=gpumatrix(gcptr,gb.dim(1),gb.dim(2),0);
                end
                return
            else
                M = ga.dim(1);
                N = ga.dim(2);
                if N ~= gb.dim(1)
                    error('matirces of incorect sizes');
                end
                K = gb.dim(2);
                if ga.mode == 0 && gb.mode ==0
                    gcptr = gpuauxSmul(ga.gptr,ga.dim(1),ga.dim(2),gb.gptr,gb.dim(1),gb.dim(2),0,0);
                    gc=gpumatrix(gcptr,M,K,0);
                else
                    gcptr = gpuauxCmul(ga.gptr,ga.dim(1),ga.dim(2),gb.gptr,gb.dim(1),gb.dim(2),0,0);
                    gc=gpumatrix(gcptr,M,K,1);
                end
            end
        end
        
        function gc = times(ga, gb)
            M = ga.dim(1);
            N = ga.dim(2);
            if M ~= gb.dim(1) || N ~= gb.dim(2)
                error('matirces of incorect sizes');
            end
            if ga.mode == 0 && gb.mode ==0
                gcptr = gpuauxSdottimes(ga.gptr,ga.dim(1),ga.dim(2),gb.gptr,gb.dim(1),gb.dim(2));
                gc=gpumatrix(gcptr,M,N,0);
            else
                gcptr = gpuauxCdottimes(ga.gptr,ga.dim(1),ga.dim(2),gb.gptr,gb.dim(1),gb.dim(2));
                gc=gpumatrix(gcptr,M,N,1);
            end

        end
        
        function gc = rdivide(ga,gb)
            M = ga.dim(1);
            N = ga.dim(2);
            if M ~= gb.dim(1) || N ~= gb.dim(2)
                error('matirces of incorect sizes');
            end
            if ga.mode == 0 && gb.mode ==0
                gcptr = gpuauxSdotdiv(ga.gptr,ga.dim(1),ga.dim(2),gb.gptr,gb.dim(1),gb.dim(2));
                gc=gpumatrix(gcptr,M,N,0);
            else
                gcptr = gpuauxCdotdiv(ga.gptr,ga.dim(1),ga.dim(2),gb.gptr,gb.dim(1),gb.dim(2));
                gc=gpumatrix(gcptr,M,N,1);
            end

            
        end
        
        function gc = horzcat(ga,gb)
            if ga.mode == 0 && gb.mode ==0
                gptr=gpuauxSconcat(ga.gptr,ga.dim(1),ga.dim(2),gb.gptr,gb.dim(2));
                gc = gpumatrix(gptr,ga.dim(1),ga.dim(2)+gb.dim(2),0);
            else
                gptr=gpuauxCconcat(ga.gptr,ga.dim(1),ga.dim(2),gb.gptr,gb.dim(2));
                gc = gpumatrix(gptr,ga.dim(1),ga.dim(2)+gb.dim(2),1);
            end
        end
        function val = subsref(ga,s)
            if strcmp(s(1).type, '()')
                if ga.mode == 0 
                    val = gpuauxSreadVal(ga.gptr,ga.dim(1),ga.dim(2), s.subs{1});
                else
                    error('gpuauxCreadVal not implemented')
                end
            elseif strcmp(s(1).type, '.')
                if strcmp(s(1).subs, 'dim')
                    val = ga.dim(s(2).subs{1});
                elseif strcmp(s(1).subs, 'mode')
                    val = ga.mode;
                elseif strcmp(s(1).subs, 'gptr')
                    val = ga.gptr;
                end
            else
                error('Illigal input to subsref');
            end
            
            
        end
        
        function gb = ctranspose(ga)
            m = ga.dim(1);
            n = ga.dim(2);

            if ga.mode == 0
                gptr=gpuauxStrans(ga.gptr,ga.dim(1),ga.dim(2));
                gb=gpumatrix(gptr,n,m,0);
            else
                gptr=gpuauxCtrans(ga.gptr,ga.dim(1),ga.dim(2));
                gb=gpumatrix(gptr,n,m,1);
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