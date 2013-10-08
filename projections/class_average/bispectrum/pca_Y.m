function [U,S,V] = pca_Y(varargin)
%PCA  Low-rank approximation in SVD form.
%
%
%   [U,S,V] = PCA(A)  constructs a nearly optimal rank-6 approximation USV'
%             to A, using 2 full iterations of the block Lanczos algorithm
%             of block size 6+2=8, started with an n x 8 random matrix,
%             when A is m x n; the reference below explains "nearly optimal."
%             The smallest dimension of A must be >= 6 when A is the only input
%             to PCA.
%
%   [U,S,V] = PCA(A,k)  constructs a nearly optimal rank-k approximation USV'
%             to A, using 2 full iterations of the block Lanczos algorithm
%             of block size k+2, started with an n x (k+2) random matrix,
%             when A is m x n; the reference below explains "nearly optimal."
%             k must be a positive integer <= the smallest dimension of A.
%
%   [U,S,V] = PCA(A,k,its)  constructs a nearly optimal rank-k approx. USV'
%             to A, using its full iterations of the block Lanczos algorithm
%             of block size k+2, started with an n x (k+2) random matrix,
%             when A is m x n; the reference below explains "nearly optimal."
%             k must be a positive integer <= the smallest dimension of A,
%             and its must be a nonnegative integer.
%
%   [U,S,V] = PCA(A,k,its,l)  constructs a nearly optimal rank-k approx. USV'
%             to A, using its full iterations of the block Lanczos algorithm
%             of block size l, started with an n x l random matrix,
%             when A is m x n; the reference below explains "nearly optimal."
%             k must be a positive integer <= the smallest dimension of A,
%             its must be a nonnegative integer,
%             and l must be a positive integer >= k.
%
%   [U,S,V] = PCA(T,Tt,params)
%   [U,S,V] = PCA(T,Tt,params,k)
%   [U,S,V] = PCA(T,Tt,params,k,its)
%   [U,S,V] = PCA(T,Tt,params,k,its,l)
%   Same as above, but takes the names of two functions that compute the
%   application of the operators T and T-transpose on a set of 
%   column vectors. "params" are any additional parameters required by T
%   and Tt. Use [] if no such paramters are required. The signature of T
%   and Tt must be either T(v) and Tt(v), or, T(params,v) and Tt(params,v).
%   When called without a column vector as its last parameter, T must
%   return m,n,cflag, where m is the number of rows in the matrix of the
%   operator, n is the number of columns, and cflag is 1 if the operator is
%   complex and 0 otherwise.  
%   
%
%   The low-rank approximation USV' is in the form of an SVD in the sense that
%   the columns of U are orthonormal, as are the columns of V, the entries
%   of S are all nonnegative, and the only nonzero entries of S appear in non-
%   increasing order on its diagonal. U is m x k, V is n x k, and S is k x k,
%   when A is m x n.
%
%   Increasing its or l improves the accuracy of the approximation USV' to A;
%   the reference below describes how the accuracy depends on its and l.
%
%
%   Note: PCA invokes RAND. To obtain repeatable results, invoke RAND('seed',j)
%         with a fixed integer j before invoking PCA.
%
%   Note: PCA currently requires the user to center and normalize the rows or
%         columns of the input matrix before invoking PCA (if such is desired).
%
%   Note: Hoping to improve its accuracy, PCA takes the adjoint of the input
%         matrix A if A has more columns than rows. However, the code will
%         produce a correct result without taking the adjoint. If taking the
%         adjoint turns out to be costly in a particular application, then the
%         user might consider deleting the lines below that take the adjoint
%         (assuming its = 0).
%
%   Note: The user may ascertain the accuracy of the approximation USV' to A
%         by invoking DIFFSNORM(A,U,S,V).
%
%
%   inputs (the first is required):
%   A -- matrix being approximated
%   k -- rank of the approximation being constructed;
%        k must be a positive integer <= the smallest dimension of A,
%        and defaults to 6
%   its -- number of full iterations of the block Lanczos algorithm to conduct;
%          its must be a nonnegative integer, and defaults to 2
%   l -- block size of the block Lanczos iterations;
%        l must be a positive integer >= k, and defaults to k+2
%
%   outputs (all three are required):
%   U -- m x k matrix in the rank-k approximation USV' to A,
%        where A is m x n; the columns of U are orthonormal
%   S -- k x k matrix in the rank-k approximation USV' to A,
%        where A is m x n; the entries of S are all nonnegative,
%        and its only nonzero entries appear in nonincreasing order
%        on the diagonal
%   V -- n x k matrix in the rank-k approximation USV' to A,
%        where A is m x n; the columns of V are orthonormal
%
%
%   Reference:
%   Vladimir Rokhlin, Arthur Szlam, and Mark Tygert,
%   A randomized algorithm for principal component analysis,
%   arXiv:0809.2274v1 [stat.CO], 2008 (available at http://arxiv.org).
%
%
%   See also PCACOV, PRINCOMP, SVDS.
%
%   Yoel Shkolnisky, October 2008
%   Based on implementation by Mark Tygert.
   

global T_func  Tt_func

if(nargin < 1)
  error('MATLAB:pca:malformedInput',...
        'There must be at least 1 input.')
end

Op=varargin{1};

if ischar(Op) % Decide if a matrix or function handle is given,
              % and parse the input accordingly.
    [T_func,Tt_func,params,k,its,l]=checkInputsT(varargin{:});
    if isempty(params) % Retrieve the dimensions of T.
        [m n,cflag] = feval(T_func); 
    else
        [m n,cflag] = feval(T_func,params); 
    end
    Aflag=0; % We are using function handles
elseif isnumeric(Op) % Matrix
    [A,k,its,l]=checkInputsA(varargin{:});    
    [m n] = size(A); % Retrieve the dimensions of A.
    Aflag=1; % We are using an explicit matrix
    
    cflag=0;
    if ~isreal(A)
        cflag=1;
    end
    
else
    error('MATLAB:pca:malformedInput',...
        'Input 1 must be a matrix or a name of a Matlab function.')
end

%
% A few more input checks
%
if(k <= 0)
  error('MATLAB:pca:malformedInput',...
        'k must be > 0.')
end

if((k > m) || (k > n))
  error('MATLAB:pca:malformedInput',...
        'k must be <= the smallest dimension of Input 1.')
end

if(its < 0)
  error('MATLAB:pca:malformedInput',...
        'its must be >= 0.')
end

if(l < k)
  error('MATLAB:pca:malformedInput',...
        'l must be >= k.')
end


%
% Check the number of outputs.
%
if(nargout ~= 3)
  error('MATLAB:pca:malformedOutput',...
        'There must be exactly 3 outputs.')
end


%
% SVD A directly if (its+1)*l >= m/1.25 or (its+1)*l >= n/1.25.
%
% if(((its+1)*l >= m/1.25) || ((its+1)*l >= n/1.25))
% 
%     if Aflag
%         if(~issparse(A))
%             [U,S,V] = svd(A,'econ');
%         end
% 
%         if(issparse(A))
%             [U,S,V] = svd(full(A),'econ');
%         end
%     else
%         if isempty(params)
%             A=feval(T_func,eye(m,n));
%         else
%             A=feval(T_func,params,eye(m,n));
%         end
%         [U,S,V] = svd(A,'econ');
%     end
% %
% % Retain only the leftmost k columns of U, the leftmost k columns of V,
% % and the uppermost leftmost k x k block of S.
% %
%   U = U(:,1:k);
%   V = V(:,1:k);
%   S = S(1:k,1:k);
% 
%   return
% 
% end


% Define functions for uniform syntax
if Aflag
    T=@(x,params) A*x;
    Tt=@(x,params) (A')*x;
    params=[];
else
    T=@T_tmp;
    Tt=@Tt_tmp;
end

%
% Take the adjoint of the matrix if m < n.
%
flag = 0;

if(m < n)
  tmpT=T;
  T=Tt;
  Tt=tmpT;

  tmpm=m;
  m=n;
  n=tmpm;

  flag = 1;
end

%
% Apply A to a random matrix, obtaining H.
%
rand('seed',rand('seed'));

if ~cflag
    H=T(2*rand(n,l)-ones(n,l),params);
else    
    H = T((2*rand(n,l)-ones(n,l)) + i*(2*rand(n,l)-ones(n,l)),params);
end

rand('twister',rand('twister'));
F = H;

%
% Apply A*A' to H a total of its times, augmenting F with the new H each time.
%
for it = 1:its
  H=Tt(H,params);
  H=T(H,params);
  F = [F H];
end

clear H;

%
% Form a matrix Q whose columns constitute an orthonormal basis
% for the columns of F.
%
[Q,R,E] = qr(F,0);
clear F R E;

%
% SVD Q'*A to obtain approximations to the singular values
% and right singular vectors of A; adjust the obtained left singular vectors
% to approximate the left singular vectors of A.
%
B=Tt(Q,params);
B=B';
[U2,S,V] = svd(B,'econ');
U = Q*U2;

clear Q U2;

%
% Retain only the leftmost k columns of U, the leftmost k columns of V,
% and the uppermost leftmost k x k block of S.
%
U = U(:,1:k);
V = V(:,1:k);
S = S(1:k,1:k);

%
% If flag is nonzero, then swap U and V.
%
if(flag ~= 0)
  flag = U;
  U = V;
  V = flag;
end

clear flag;

end

%% checkInputsT
% First argument is a function handle. Functions for computing the operator
% and its transpose are given. Parse the input arguments accordingly.
function [T_func,Tt_func,params,k,its,l]=checkInputsT(varargin)

narg=numel(varargin);

% Check the number of inputs.

if(narg < 3)
    error('MATLAB:pca:malformedInput',...
        'When providing functions, there must be at least 3 inputs. Use [] if needed.')
end

if(narg > 6)
    error('MATLAB:pca:malformedInput',...
        'When providing functions, there must be at most 6 inputs.')
end

T_func=varargin{1};
Tt_func=varargin{2};
params=varargin{3};

% Check arguments type
if ~ischar(T_func) || ~ischar(Tt_func)
    error('MATLAB:pca:malformedInput',...
        'First two input must be the names of the operator and its adjoint.')
end
    
% Set the inputs k, its, and l to default values, if necessary.
if(narg == 3)
    k=6;
else
    k=varargin{4};
    assertInt(k,4);
end

if(narg <= 4)
    its = 2;
else
    its=varargin{5};
    assertInt(its,5);
end

if(nargin <= 5)
    l = k+2;
else
    l=varargin{6};
    assertInt(l,6);
end

end

%% checkInputsA
% The matrix A is given explicitly. Parse the input arguments accordingly.
function [A,k,its,l]=checkInputsA(varargin)

narg=numel(varargin);

% Check the number of inputs.
if(narg > 4)
    error('MATLAB:pca:malformedInput',...
        'When providing a matrix, there must be at most 4 inputs.')
end

A=varargin{1};

% Check argument type
if(~isnumeric(A))
    error('MATLAB:pca:malformedInput',...
        'Input 1 must be a floating-point matrix.')
end

if(isempty(A))
    error('MATLAB:pca:malformedInput',...
        'Input 1 must not be empty.')
end

% Set the inputs k, its, and l to default values, if necessary.
if(narg == 1)
    k=6;
else
    k=varargin{2};
    assertInt(k,2);
end

if(narg <= 2)
    its = 2;
else
    its=varargin{3};
    assertInt(its,3);
end

if(nargin <= 3)
    l = k+2;
else
    l=varargin{4};
    assertInt(l,4);
end


end

%% assertInt
% Check that the argument is an integer scalar
function assertInt(x,pos)
if ~isscalar(x) || ~isreal(x) || (floor(x)~=x)
    error('MATLAB:pca:malformedInput',...
        'Input %d must be an integer',pos);
end
end

%% T_tmp
function w=T_tmp(v,params)
global T_func
if isempty(params)
    w=feval(T_func,v);
else
    w=feval(T_func,params,v);
end
end

%% Tt_tmp
function w=Tt_tmp(v,params)
global Tt_func
if isempty(params)
    w=feval(Tt_func,v);
else
    w=feval(Tt_func,params,v);
end
end