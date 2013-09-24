%
%  [x,y,info,r] = sdplr(A,b,c,K,pars,lrA,x0,y0,info0,r0)
%
%  SDPLR 1.03-beta   http://dollar.biz.uiowa.edu/~sburer/software/SDPLR
%                    Email and bug reports to: samuel-burer@uiowa.edu
%
%  > x = sdplr(A,b,c,K)
%
%    Solves the semidefinite program:   min  c'*x  st  A*x = b, x in K
%
%    K   describes  the   cone  of   nonnegative  linear   and  positive
%    semidefinite variables. K follows the SeDuMi format (see the SeDuMi
%    website http://sedumi.mcmaster.ca/).
%
%    K  can  have two  fields,  K.l  and  K.s, representing  Linear  and
%    Semidefinite. K.l is the number of nonnegative components of x, and
%    K.s lists the dimensions of the semidefinite constraints on x.
%
%    Example: K.l = 10,  K.s = [4 3], and K.l listed  before K.s. Then x
%    has dimension 35 = 10 + 4*4 + 3*3 and satisfies the constraints
%
%                   x(1:10)                nonnegative
%                   reshape(x(11:26),4,4)  positive semidefinite
%                   reshape(x(27:35),3,3)  positive semidefinite
%
%    Example: K.s = [3  4], K.l = 10, and K.s listed  before K.l. Then x
%    has dimension 35 = 10 + 4*4 + 3*3 and satisfies the constraints
%
%                   reshape(x(1:9),3,3)    positive semidefinite
%                   reshape(x(10:25),4,4)  positive semidefinite
%                   x(26:35)               nonnegative
%
%  > [x,y] = sdplr(A,b,c,K)
%
%    Returns the dual solution of:   max  b'*y  st  c - A'*y in K
%
%  > [x,y,info] = sdplr(A,b,c,K)
%
%    Returns information from the execution of SDPLR, where:
%
%        info.majoriter  =  # of major iterations (i.e., penalty increases)
%        info.minoriter  =  # of minor iterations (i.e., steps)
%        info.cgs        =  # of conjugate gradient iterations
%        info.time       =  total time (in seconds)
%        info.penalty    =  final value of penalty parameter
%
%  > [x,y,info,r] = sdplr(A,b,c,K)
%
%    Returns the internal format r of x used by SDPLR.
%
%    The variable r  is a cell array, such that  r{i} stores the portion
%    of x correpsonding to the i-th cone constraint on x, as given by K:
%
%      portion x  =  r{i}.*r{i}                    if linear
%      portion x  =  reshape(r{i}*r{i}',sz*sz,1)   if semidef of size sz
%
%    Example: K.l = 10, K.s = [4 3], and K.l listed before K.s. Then
%
%                   x(1:10)   =  r{1}.*r{1}
%                   x(11:26)  =  reshape(r{2}*r{2}',16,1)
%                   x(27:35)  =  reshape(r{3}*r{3}',9,1)
%
%  > [x,y,info,r] = sdplr(A,b,c,K,pars)
%
%    Passes parameters into SDPLR, where:
%
%        pars.feastol    = Desired level of (scaled) infeasibility in
%                          A*x = b (default = 1.0e-5)
%
%        pars.centol     = Desired level of accuracy in intermediate
%                          calculations (default = 1.0e-1)
%
%        pars.dir        = Method for calculating step direction
%                          (1 = limited memory BFGS, 2 = truncated Newton)
%                          (default = 1)
%
%        pars.penfac     = Factor by which the penalty parameter is
%                          increased each major iteration (default = 2.0)
%
%        pars.reduce     = Whether or not to reduce the problem dimension-
%                          ality as the algorithm progresses (1 = yes,
%                          0 = no) (default = 0)
%
%        pars.limit      = Limit on the total time (in seconds)
%                          (default = 3600)
%
%        pars.printlevel = Whether to print output (0 = no, 1 = yes)
%                          (default = 1)
%
%        pars.seed       = Seed for random number generator.
%
%        pars.forcerank  = Array having  the same  number of  entries as
%                          the total  number of  blocks  specified by  K
%           (where each diagonal block counts  as one block). Each array
%           entry corresponds to a block of x (in order) and tells SDPLR
%           to force a rank condition on that block. Diagonal blocks can
%           only have rank  1, while semidefinite blocks of  size sz can
%           have any rank <= sz.
% 
%           Example: K.l =  10, K.s = [4 3], and  K.l listed before K.s.
%           Then the array [1 2 2] forces the two semidefinite blocks to
%           be rank 2.
%
%           Warning:  When   the  rank   is  restricted,   SDPLR  cannot
%           theoretically guarantee  convergence to an  optimal solution
%           of the rank-restricted problem!
%
%   To save on memory, there is also:
%
%        pars.soln_factored  =  Return x in its factored form (i.e., x <- r)
%                               (1 = yes, 0 = no) (default = 0)
%
%  > [x,y,info,r] = sdplr(A,b,c,K,pars,lrA)
%
%    Passes low-rank data matrices into SDPLR,  where lrA is an array of
%    structures.  E.g.,  lrA =  [lrA(1),  lrA(2),  lrA(3)] passes  three
%    low-rank matices.
%
%    Each low-rank  matrix must apply to  a (semidefinite) cone i  and a
%    constraint j,  where the cones are  given by K and  the constraints
%    correspond to  the rows of A.  One can also indicate  the objective
%    function by j=0. The matrix is  given by V*diag(D)*V', where V is a
%    matrix and D  is a vector. The number  of rows of V is  the size of
%    the semidefinite cone, and the number of columns of V (also size of
%    D) is the rank of the matrix.
%
%    Each structure in lrA has four fields:
%
%        cons   =  which constraint j
%        start  =  index of the first position of cone i in x
%        D      =  vector D
%        V      =  matrix V
%
%    Example: K.l = 10, K.s = [4 3],  and K.l listed before K.s. If V is
%    size (3,2) and D is size 2, then the rank-2 matrix V*diag(D)*V' can
%    be applied to  x in the third cone (semidef)  and second constraint
%    as
%
%         lrA(1).cons   =  2 
%         lrA(1).start  =  27
%         lrA(1).D      =  D
%         lrA(1).V      =  V
%
%    Note: Even if  the j-th constraint is specified  completely by lrA,
%    it is still necessary  for A to have a j-th  row (sparse and empty)
%    to serve as a placeholder.
%
%  > [x,y,info,r] = sdplr(A,b,c,K,pars,lrA,x0,y0,info0,r0)
%
%    Passes   initial  state   of  algorithm   into  SDPLR.   Format  of
%    (x0,y0,info0,r0) is exactly that of output (x,y,info,r). Any subset
%    of initial information may be specified. However, if both x0 and r0
%    are specified  then SDPLR reads  only r0. Also,  x0 or r0  may need
%    to  be  adjusted to  match  the  internal rank-restricted  positive
%    semidefinite variable of SDPLR.
%
%  ------------------------------------------------------------------------
%
%  SDPLR 1.03-beta   http://dollar.biz.uiowa.edu/~sburer/software/SDPLR
%                    Email and bug reports to: samuel-burer@uiowa.edu
function [varargout] = sdplr(varargin);

% If initial information has been given...
if length(varargin) >= 7
  % First allocate recognizable names
  b = varargin{2};
  K = varargin{4};
  pars = varargin{5};
  x0 = varargin{7};
  if length(varargin) >=  8; y0    = varargin{ 8}; else; y0    = []; end;
  if length(varargin) >=  9; info0 = varargin{ 9}; else; info0 = []; end;
  if length(varargin) >= 10; r0    = varargin{10}; else; r0    = []; end;
  % Setup some structures that may be used below (also do a bit of error checking)
  fields = fieldnames(K);
  if length(fields) > 2
    fprintf('Error: K has too many fields.\n');
    return
  end
  for i = 1:length(fields)
    if fields{i} ~= 's' & fields{i} ~= 'l'
      fprintf('Error: K has illegal field name.\n');
      return
    end
  end
  numblk = 0;
  if isfield(K,'l')
    numblk = numblk + 1;
  end
  if isfield(K,'s')
    numblk = numblk + length(K.s);
  end
  if max(size(y0)) > 0 & max(size(y0)) ~= max(size(b))
    fprintf('Error: size of y0 does not match size of b.\n');
    return
  end
  % If both x0 and r0 have been given, then print warning but drop x0
  if max(size(x0)) > 0 & max(size(r0)) > 0
    fprintf('Warning: Both inputs 7 and 10 have been given. Ignoring input 7.\n');
    x0 = [];
  end
  % If x0 has been given (and r0 has not), then convert x0 to r0 and drop x0
  if max(size(x0)) > 0 & max(size(r0)) == 0
    % Extract field names; will allow us to determine which of K.l or
    % K.s comes first. (Also perform some error checking)
    pos = 1; blk = 1;
    if fields{1} == 'l'
      r0{blk} = mysqrt( x0( pos : pos + K.l - 1 ) , 'vec' ); blk = blk+1; pos = pos + K.l;
      if length(fields) > 1
        for j = 1:length(K.s)
          r0{blk} = mysqrt( x0( pos : pos + (K.s(j))^2 - 1 ) , 'mat' ); blk = blk+1; pos = pos + (K.s(j))^2;
        end
      end
    end
    if fields{1} == 's'
      for j = 1:length(K.s)
        r0{blk} = mysqrt( x0( pos : pos + (K.s(j))^2 - 1 ) , 'mat' ); blk = blk+1; pos = pos + (K.s(j))^2;
      end
      if length(fields) > 1
        r0{blk} = mysqrt( x0( pos : pos + K.l - 1 ) , 'vec' ); blk = blk+1; pos = pos + K.l;
      end
    end
    x0 = [];
  end
  % If pars.forcerank is available, then go ahead and chop r0 if
  % necessary (extending r0 will be done inside C portion of code)
  if isfield(pars,'forcerank')
    if max(size(pars.forcerank)) == numblk & min(size(pars.forcerank)) == 1
      blk = 1;
      if fields{1} == 'l'
        blk = blk+1;
        if length(fields) > 1
          for j = 1:length(K.s)
            if size(r0{blk},2) > pars.forcerank(blk)
              r0{blk} = r0{blk}(:,1:pars.forcerank(blk));
            end
            blk = blk+1;
          end
        end
      end
      if fields{1} == 's'
        for j = 1:length(K.s)
          if size(r0{blk},2) > pars.forcerank(blk)
            r0{blk} = r0{blk}(:,1:pars.forcerank(blk));
          end
          blk = blk+1;
        end
        if length(fields) > 1
          blk = blk+1;
        end
      end
    else
      fprintf('Warning: forcerank contains incorrect or incomplete information.\n');
    end
  end
end

% Put arguments back in place (perhaps more than original)
if length(varargin) >= 7
  varargin{ 7} = x0;
  varargin{ 8} = y0;   
  varargin{ 9} = info0;
  varargin{10} = r0;
end

if nargout == 0
    mexsdplr(varargin{:});
else
    [varargout{1:nargout}] = mexsdplr(varargin{:});
end

return

function result = mysqrt(input, type)

if type == 'vec'
  result = sqrt(max(input,0));
end
if type == 'mat'
  dim = sqrt(max(size(input)));
  mat = reshape(input,dim,dim);
  mat = 0.5*(mat + mat');
  [V,D] = eig(mat);
  d = sqrt(max(diag(D),0));
  ind = find(d > 0);
  d = d(ind);
  V = V(:,ind);
  [d,ind] = sort(d,1,'descend');
  V = V(:,ind);
  result = V*diag(d);
end
return
