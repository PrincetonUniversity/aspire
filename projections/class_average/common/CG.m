% function [Y,flag,relres,iter,absres] = CG(PtP,X,params,ErrTol,MaxIts,guess,verbose,RefY)
%
% Solve the system X=PtP(Y,params) using the conjugate gradient method.
% The operator PtP must be hermitian and positive defined.
% The firat parameter to PtP must be the variable Y. params can be any list
% of comma separated arguments.
%
%  Input parameters:
%    PtP       Name of the operator to invert
%    X         The transformed matrix. The matrix at the range space of the operator PtP, 
%              whose source needs to be found.
%    params    Additional parameters to the operator PtP.
%    ErrTol    Error tolerance of the CG method. Default 1.e-9.
%    MaxIts    Maximum number of iterations. Default 10.
%    guess     Initial guess of the solution. Default is X.
%    versos    By default, if more than one output argument is specified, then all output 
%              messages are suppressed. Set this flag to any value other than 0 to
%              always display output messages.
%    RefY      The untransformed matrix Y. Used only for checking absolute error.
%              If not specified, absolute error is not computed.
%
%  Output parameters:
%    Y         The result matrix of the CG method. This is the estimate of the CG
%              method to the solution of the system X=PtP(Y).
%    flag      A flag that describes the convergence of the CG method.
%              0 CG converged to the desired tolerance ErrTol within MaxIts iterations.
%              1 CG did not converge to ErrTol within MaxIts iterations.
%    relres    Residual error at the end of the CG method. Computed using max norm.
%    iter      The iteration number at which ErrTol was achieved. Relevant only if
%              flag=0.
%    absres    The absolute error at the end of the CG method. Relevant only if RefY
%              was given as parameter. Computed using max norm.
%
% Yoel Shkolnisky 15/12/02

function [Y,flag,relres,iter,absres] = CG(PtP,X,params,ErrTol,MaxIts,guess,verbose,RefY)

% Check the input and initialize flags and default parameters
ref_given=1;         % The flags is 1 if the reference untransformed matrix RefY is given and 0 otherwise

if nargin<8         % RefY not given
   ref_given=0;   
end

if nargin<7
   verbose=0;
end

if nargin<6        % No initisl guess given
    guess=0;
end
    
if nargin<5        % MaxIts not given      
   MaxIts = 10;
end

if nargin<4        % ErrTol not given
    ErrTol = 1.e-9;
end

if nargin<3        % no additional parameters to PtP are given
    params=cell(0);
end

% Initialize convergence flag. If the routine will detect that the CG method converged, this flag 
% will be set to 0 to represent convergence. By default it assumes that the CG did not converge.
flag=1;

% Set flag to suppress output if "flag" output is specified.
suppress_output=0;
if (nargout>1) & (~verbose)
   suppress_output=1;
end

% iter holds the iteration in which CG converged to ErrTol.
iter=-1;

% Initialization
xk = guess;
temp = feval(PtP,xk,params{:});
gk = temp-X;
pk = -gk;
dk = sum(abs(gk(:)).^2);

% Conjugate gradient iteration
j=2;
done = 0;
while (j<=MaxIts) & (~done)
    perr=sum((abs(pk(:))).^2);
    if ref_given % If reference matrix is given compute absolute error
        xerr=max(max(abs(RefY-xk)));
    end
    if (~suppress_output) & (flag)
        fprintf('Iteration %2d:  Residual error=%-2.7e',j-1,perr);
        if ref_given
            fprintf('\t Absolute error=%-2.7e',xerr);
        end
        fprintf('\n');
    end

    if perr<=ErrTol
        iter=j-1;  %CG converged at previous iteration
        flag=0;
        done=1;
    end
    if perr>ErrTol
        hk = feval(PtP,pk,params{:});
        tk = dk/dot(pk(:),hk(:));  %line search parameter
        xk = xk+tk*pk;       % update approximate solution
        gk = gk+tk*hk;       % update gradient
        temp = sum(abs(gk(:)).^2);
        bk = temp/dk;
        dk = temp;
        pk = -gk+bk*pk;       %update search direction
    end
    j=j+1;
end

relres = perr;
if ref_given
    absres = xerr;
end      

Y = xk;
