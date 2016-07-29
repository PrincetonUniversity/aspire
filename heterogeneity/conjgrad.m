% CONJGRAD Solve system using conjugate gradients
%
% Usage
%    [x, obj, info] = conjgrad(b, Afun, cg_opt, init);
%
% Input
%    b: The vector consisting of the right hand side of Ax = b.
%    Afun: A function handle specifying the linear operation x -> Ax.
%    cg_opt: The parameters for the conjugate gradient method, including:
%       max_iter: Maximum number of iterations (default 10).
%       verbose: The extent to which information on progress should be output
%          to the terminal (default 1).
%       iter_callback: If non-empty, specifies a function to be called at the
%          end of every iteration. In this case, iter_callback must be a
%          function handle taking as single argument the info structure at
%          the current iteration. For information on the info structure, see
%          below (default []).
%       preconditioner: If non-empty, specifies a preconditioner to be used in
%          every iteration as a function handle defining the linear operator
%          x -> Px (default []).
%       rel_tolerance: The relative error at which to stop the algorithm, even
%          if it has not yet reached the maximum number of iterations (default
%          1e-15).
%       store_iterates: Defines whether to store each intermediate results in
%          the info structure under the x field. Since this may require a
%          large amount of memory, it is recommended not to do this unless
%          necessary (default false).
%    init: A structure specifying the starting point of the algorithm. This
%       can contain values of x or p that will be used for intialization
%       (default empty).
%
% Output
%    x: The result of the conjugate gradient method after max_iter iterations
%       or once the residual norm has decreased below rel_tolerance, relative.
%    obj: The value of the objective function at the last iteration.
%    info: A structure array containing intermediate information obtained
%       during each iteration. These fields include:
%       iter: The iteration number.
%       x (if store_iterates is true): The value of x.
%       r: The residual vector.
%       p: The p vector.
%       res: The square norm of the residual.
%       obj: The objective function.

function [x, obj, info] = conjgrad(b, Afun, cg_opt, init)
    if nargin < 3
        cg_opt = struct();
    end

    if nargin < 4
        init = struct();
    end

    if ~isfield(init, 'x') || isempty(init.x)
        x = zeros(size(b));
    else
        x = init.x;
    end

    if ~isfield(cg_opt, 'max_iter'), cg_opt.max_iter = 10; end
    if ~isfield(cg_opt, 'verbose'), cg_opt.verbose = 1; end
    if ~isfield(cg_opt, 'iter_callback'), cg_opt.iter_callback = []; end
    if ~isfield(cg_opt, 'preconditioner'), cg_opt.preconditioner = @(x)(x); end
    if ~isfield(cg_opt, 'rel_tolerance'), cg_opt.rel_tolerance = 1e-15; end
    if ~isfield(cg_opt, 'store_iterates'), cg_opt.store_iterates = false; end

    b_norm = sqrt(sum(abs(b).^2, 1));
    
    r = b;
    s = cg_opt.preconditioner(r);
    if any(x(:) ~= 0)
        if cg_opt.verbose
            fprintf('[CG] Calculating initial residual...');
        end
        ticker = tic;
        Ax = Afun(x);
        r = r-Ax;
        s = cg_opt.preconditioner(r);
        if cg_opt.verbose
            fprintf('OK (%.2f s)\n', toc(ticker));
        end
    else
        Ax = zeros(size(x));
    end

    old_res = real(sum(conj(r).*s, 1));

    obj = real(sum(conj(x).*Ax, 1)) - 2*real(sum(conj(b).*x, 1));

    if ~isfield(init, 'p') || isempty(init.p)
        p = s;
    else
        p = init.p;
    end

    info = struct();

    info(1).iter = 0;
    if cg_opt.store_iterates
        info(1).x = x;
    end
    info(1).r = r;
    info(1).p = p;
    info(1).res = old_res;
    info(1).obj = obj;

    if cg_opt.verbose
        fprintf(['[CG] Initialized. Residual: %g. ' ...
            'Objective: %g.\n'], sum(info(1).res), sum(info(1).obj));
    end

    for iter = 1:cg_opt.max_iter
        if cg_opt.verbose
            fprintf('[CG] Applying matrix & preconditioner...');
        end
        ticker = tic;
        Ap = Afun(p);

        alpha = old_res./real(sum(conj(p).*Ap, 1));
        x = x + bsxfun(@times, alpha, p);
        Ax = Ax + bsxfun(@times, alpha, Ap);

        r = r - bsxfun(@times, alpha, Ap);
        s = cg_opt.preconditioner(r);
        new_res = real(sum(conj(r).*s, 1));
        beta = new_res./old_res;
        p = s + bsxfun(@times, beta, p);

        if cg_opt.verbose
            fprintf('OK (%.2f s)\n', toc(ticker));
        end
        
        obj = real(sum(conj(x).*Ax, 1)) - 2*real(sum(conj(b).*x, 1));
    
        old_res = new_res;

        info(iter+1).iter = iter;
        if cg_opt.store_iterates
            info(iter+1).x = x;
        end
        info(iter+1).r = r;
        info(iter+1).p = p;
        info(iter+1).res = old_res;
        info(iter+1).obj = obj;

        if cg_opt.verbose
            fprintf(['[CG] Iteration %d. Residual: %g. ' ...
                'Objective: %g.\n'], iter, sum(info(iter+1).res), sum(info(iter+1).obj));
        end
    
        if ~isempty(cg_opt.iter_callback)
            cg_opt.iter_callback(info);
        end

        if all(sqrt(new_res) < b_norm*cg_opt.rel_tolerance)
            break;
        end
    end
end
