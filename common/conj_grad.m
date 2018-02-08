% CONJ_GRAD Solve system using conjugate gradients
%
% Usage
%    [x, obj, info] = conj_grad(Afun, b, cg_opt, init);
%
% Input
%    Afun: A function handle specifying the linear operation x -> Ax.
%    b: The vector consisting of the right hand side of Ax = b.
%    cg_opt: The parameters for the conjugate gradient method, including:
%          - max_iter: Maximum number of iterations (default 50).
%          - verbose: The extent to which information on progress should be
%             output to the terminal (default 1).
%          - iter_callback: If non-empty, specifies a function to be called at
%             the end of every iteration. In this case, iter_callback must be a
%             function handle taking as single argument the info structure at
%             the current iteration. For information on the info structure, see
%             below (default []).
%          - preconditioner: If non-empty, specifies a preconditioner to be
%             used in every iteration as a function handle defining the linear
%             operator x -> Px (default []).
%          - rel_tolerance: The relative error at which to stop the algorithm,
%             even if it has not yet reached the maximum number of iterations
%             (default 1e-15).
%          - store_iterates: Defines whether to store each intermediate results
%             in the info structure under the x, p and r fields. Since this
%             may require a large amount of memory, this is not recommended
%             (default false).
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
%          - iter: The iteration number.
%          - x (for store_iterates true): The value of x.
%          - r (for store_iterates true): The residual vector.
%          - p (for store_iterates true): The p vector.
%          - res: The square norm of the residual.
%          - obj: The objective function.

function [x, obj, info] = conj_grad(Afun, b, cg_opt, init)
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

    cg_opt = fill_struct(cg_opt, ...
        'max_iter', 50, ...
        'verbose', 0, ...
        'iter_callback', [], ...
        'preconditioner', @(x)(x), ...
        'rel_tolerance', 1e-15, ...
        'store_iterates', false);

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

    old_res = sqrt(real(sum(conj(r).*s, 1)));

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
        info(1).r = r;
        info(1).p = p;
    end
    info(1).res = sqrt(sum(abs(r).^2, 1));
    info(1).obj = obj;

    if cg_opt.verbose
        fprintf(['[CG] Initialized. Residual: %g. ' ...
            'Objective: %g.\n'], norm(info(1).res), sum(info(1).obj));
    end

    if b_norm == 0
        return;
    end

    iter = 1;
    while true
        if iter >= cg_opt.max_iter
            break;
        end

        if cg_opt.verbose
            fprintf('[CG] Applying matrix & preconditioner...');
        end
        ticker = tic;
        Ap = Afun(p);

        old_gamma = real(sum(conj(s).*r, 1));

        alpha = old_gamma./real(sum(conj(p).*Ap, 1));
        x = x + bsxfun(@times, alpha, p);
        Ax = Ax + bsxfun(@times, alpha, Ap);

        r = r - bsxfun(@times, alpha, Ap);
        s = cg_opt.preconditioner(r);
        new_gamma = real(sum(conj(r).*s, 1));
        beta = new_gamma./old_gamma;
        p = s + bsxfun(@times, beta, p);

        if cg_opt.verbose
            fprintf('OK (%.2f s)\n', toc(ticker));
        end

        obj = real(sum(conj(x).*Ax, 1)) - 2*real(sum(conj(b).*x, 1));

        old_gamma = new_gamma;

        info(iter+1).iter = iter;
        if cg_opt.store_iterates
            info(iter+1).x = x;
            info(iter+1).r = r;
            info(iter+1).p = p;
        end
        info(iter+1).res = sqrt(sum(abs(r).^2, 1));
        info(iter+1).obj = obj;

        if cg_opt.verbose
            fprintf(['[CG] Iteration %d. Residual: %g. Objective: %g.\n'], ...
                iter, norm(info(iter+1).res), sum(info(iter+1).obj));
        end

        if ~isempty(cg_opt.iter_callback)
            cg_opt.iter_callback(info);
        end

        if all(info(iter+1).res < b_norm*cg_opt.rel_tolerance)
            break;
        end

        iter = iter+1;
    end

    if iter == cg_opt.max_iter
        warning('Conjugate gradient reached maximum number of iterations!');
    end
end
