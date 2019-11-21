% compute the spherical bessel coefficients of a cartesian volume using
% gradient descent for the cost function, with step size found by exact
% line search.

% input:
% V: cartesian volume
% L: maximum bandlimit

% output: 
% x: vector of coefficients

% Yuan Liu, 7/17/2017

% cd /u/liuyuan/Documents/MATLAB/Subtomo_average/ASPIRE
% initpath;
% close all;
% format long;
R = 32;
L = R;

%% Generate synthetic example: Non-centered Gaussian from Boris
x_1d = (-L:1:L)/L;   % - Odd number of points
[x_3d,y_3d,z_3d] = meshgrid(x_1d,x_1d,x_1d);
r_b = sqrt(x_3d.^2 + y_3d.^2 + z_3d.^2);
ball = (r_b <= 1);

t = 0.05;
vol = exp(-((x_3d-0.3).^2 + (y_3d-0.2).^2 + z_3d.^2)/t);
hatv = vol;

%% 70S E.coli ribosome
% load('cleanrib.mat')
% 
% maxL = 20;%31;
% hatv = volref(1:end-1,1:end-1,1:end-1);
r_mask = 2;
r_select_ratio = R;

MACHINEEPS = 1e-15;
mask_func = @(t) 1;
abs_flag = 0;
mapping_f= @(z) z;

% makd 3d grid
N = (size(hatv,1)-1)/2;
[x, y, z]=meshgrid(-N:N, -N:N, -N:N);
x = x/N;
y = y/N;
z = z/N;
% convert to r, theta, phi
r = sqrt( x.^2 + y.^2 + z.^2 );
theta = acos( z./r);
phi = acos( x./(r.*sin(theta)) );
j1 = find( y < 0);
phi(j1) = 2*pi - phi(j1) ;
j = find( theta < MACHINEEPS | theta > pi-MACHINEEPS); % z = 1 or -1, x=y=0
phi(j) = 0;
phi = real( phi); %enforce real numbers
jorigin = find( r < MACHINEEPS ); % be careful at orgin point
theta(jorigin) = 0;
phi( jorigin ) = 0;

jball = find( r < r_mask ); % semiball, compact support on r < r_mask
jjorigin = find( r(jball) < MACHINEEPS );

hatv = mapping_f( hatv ); %pointwise filter
radial_mask = mask_func( r);
if abs_flag
    hatv_m = abs(hatv).*radial_mask; % delete it???
else
    hatv_m = hatv.*radial_mask;
end
y = hatv_m( jball);

figure(4);
plot(y);

%% load pre-computed roots of spherical bessel
tic
load SphericalBessel.mat
B = bessel;
%r_mask = N;
B = B( B(:, 3)<pi * r_mask * r_select_ratio & B(:,1) <= maxL, :); %Nyquist criterion
l_grid = B(:, 1);
s_grid = B(:, 2);
R_ls   = B(:, 3); %zeros of j_l(r)
clear B bessel;

% make the basis function
siz_grid = numel( jball);

%%% basis of Psi_lms
jl_func = @(l, x) besselj(l+1/2, x).*sqrt(pi./(2*x)); % x > 0

num_basis = 0;
for ll = 0: maxL
    ind_ls = find( l_grid == ll);
    siz_ls = numel( ind_ls );
    num_basis = num_basis + siz_ls*(2*ll+1);
end
tablelms = zeros( 3, num_basis); % the three rows are l, m, s
%Psilms = zeros( siz_grid, num_basis); %store the basis
indcol = 0;

fprintf('generating basis:');
for ll = 0:maxL   
    fprintf('%d, ', ll);
     
    %%% spherical harmonics for order ll
    [X, ~, dems]= xlm( ll, -ll:1:ll, theta(jball));    
    X = X.'; %each column is a basis
    Y = exp( sqrt(-1)* phi(jball)* (-ll:1:ll));
    Y = Y.*X; % Y_l^m( theta, phi) = X_lm( theta )*exp(i*m*phi)
              % Y_l^{-m} = (-1)^m * conj( Y_l^m ),
              % X_{l,-m}(theta) = (-1)^m*X_{lm}(theta)
    
    % be careful with origin point, later in the radial part only when ll =0 it can be a constant at the origin
    if ll == 0 
        Y(jjorigin ) = Y(1, 1); %Y_00 is the constant
    else       
        Y(jjorigin, :) = 0;
    end
%     if sum(isnan(Y(:)) + isinf(Y(:))) > 0
%         error('ylm take value NaN or Inf.')
%     end
    
    %%% spherical bessel series, Jls = j_l(R_ls*r)
    ind_ls = find( l_grid == ll);
    siz_ls = numel( ind_ls );
    Jls = zeros( siz_grid , siz_ls);
    for ss = 1: siz_ls
        root_ls = R_ls( ind_ls( ss )); % s-th root of j_l(r)
        normalcont_ls = sqrt(2)/abs( jl_func(ll+1, root_ls) );%/sqrt(pi^3);
        
        Jls(:,ss) = jl_func( ll, root_ls*( r(jball)/r_mask ))*normalcont_ls;    
        % be careful with the origin point
        if ll == 0
            Jls(jjorigin, ss) = normalcont_ls;
        else
            Jls(jjorigin, ss) = 0;
        end       
    end
%     if sum(isnan(Jls(:)) + isinf(Jls(:))) > 0
%         error('Jls take value NaN or Inf.')
%     end
   
    for ii = 1: (2*ll+1)
        tablelms(1, indcol+1:indcol+siz_ls) = ll; %l
        tablelms(2, indcol+1:indcol+siz_ls) = dems(ii);%m
        tablelms(3, indcol+1:indcol+siz_ls) = 1:siz_ls; %s        
        %Psilms(  :, indcol+1:indcol+siz_ls) = bsxfun(@times, Jls, Y(:,ii)); 
                                              %orthonormal basis in the unit ball        
        indcol = indcol + siz_ls;
    end
end
fprintf('\n');
toc

%% find coefficients with gradient descent

% hatv = zeros(siz_grid,1);
% y = hatv;
% termination tolerance
tol = MACHINEEPS*10;

% maximum number of allowed iterations
maxiter = 50;
D = ones(maxiter,1);
%D = cat(1,D,zeros(maxiter,1));

% minimum allowed perturbation
dxmin = 1e-14;

% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; 
x0 = 1e-4*ones(num_basis,1);
%x0 = a;
x = x0; 
niter = 0; 
dx = inf;
dD = 1;

% gradient descent algorithm:
tic
while and(gnorm>=tol, and(niter+1 <= maxiter, dD >= tol))
    
    % calculate gradient:
    g = 2*Psilms'*(Psilms*x-y);
    gnorm = norm(g);
    
    % record status
    D(niter+1) = norm(Psilms*x-y)^2;
    if niter == 0
        dD = D(1);
    else 
        dD = abs(D(niter+1) - D(niter));
    end
    
    fprintf('niter:%d, dx:%.14f, gnorm:%.14f, dobjective:%.14f\n', niter, dx, gnorm, dD);
    
    % take step:
    t = (Psilms'*(Psilms*x-y))./(Psilms'*Psilms*g);
    xnew = x - t.*g;
    
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    
    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;    
end
t_gd = toc
xopt = x;
%fopt = f2(xopt);
niter = niter - 1;

%%
figure(3)
semilogy(1:niter+1,D(1:niter+1))
xlabel('iteration','FontSize',16);
ylabel('objective value','FontSize',16);
title('Gradient descent iterations for clean Ribosome 80S volume','FontSize',16);

%% get clm from x
clm = cell( maxL+1 , 1);

for ll = 0 : maxL
    
    ind_ls = find( l_grid == ll);
    siz_ls = numel( ind_ls );   
    al = reshape( x(  tablelms(1,:) == ll ), siz_ls, 2*ll+1);
    
    % A_{l,-m} = (-1)^(l+m)* conj( A_{lm} ), thus A_l0 = (-1)^l*conj( A_l0)
    if mod(ll, 2) == 0
        al(:, ll+1) = real( al(:, ll+1) );
    else
        al(:, ll+1) = sqrt(-1)*imag( al(:, ll+1) );
    end
    
    clm{ ll +1 } = al;
end

   
%% plot Alm(r)
figure(1), clf,
%xx = 0.01: r_mask/100 :r_mask;
%xx = (0:0.1:N)/N;
xx = (0:N-1)/N;
f = zeros( 1, numel( xx));

for ll = 0: 1 : maxL
   
    [X, ~, dems]= xlm( ll, -ll:1:ll, pi/2);
    X = X';
    Y = exp( sqrt(-1)* pi/2* (-ll:1:ll));
    Y = Y.*X;
      
   for mm = -ll:ll  %only plot nonnegative m, A_{l,-m} = (-1)^(l+m)* conj( A_{lm} )
            
       a_llmm = clm{ ll+1}(:, mm+ll+1);
       
       ind_ls = find( l_grid == ll);
       siz_ls = numel( ind_ls);
       if size( a_llmm, 1) ~= siz_ls
           error('size of a_llmm not equal to siz_ls');
       end
       
       for ss = 1: siz_ls
           root_ls = R_ls( ind_ls( ss )); % s-th root of j_l(r)
           normalcont_ls = sqrt(2)/abs( jl_func(ll+1, root_ls) );
           jls = jl_func( ll, root_ls*(xx/r_mask))*normalcont_ls;
           jls = jls*Y(mm+ll+1);
           if ll == 0
               jls(1) = normalcont_ls*Y(1);
           else
               jls(1) = 0;
           end
           f = f + a_llmm(ss)*jls;
       end
       
   end
end

%F_x = Psilms(jjorigin:jjorigin+N-1,:)*x;
%F_r = mvnpdf([(0:0.1:N)' zeros(10*N+1,1) zeros(10*N+1,1)],mu1,Sigma1);
set(gca, 'FontSize', 16);
hold on;
%plot( xx, real(f),'o-',xx, imag(f), '--g', xx, F_r','+-');
plot( xx, real(f),'o-',xx, imag(f), '--g', xx, y(jjorigin:jjorigin+N-1),'+-');
title(sprintf('Gradient descent interpolation for radial with maxL=%d', ll));
legend('approximation (real part)','approximation (complex part)','original function')
xlabel('radius R/Rmask');
ylabel('function value');

%% least-square
fprintf('(|lms|, |omega_grid|) = (%d, %d), \nsolving ls...\n', num_basis, siz_grid);

%b = hatv(jball);
b = hatv_m( jball);
disp('Computing A_lm with least squares')
tic
a = (Psilms'*Psilms)\(Psilms'*b);       %a = Psilms\b;
t_ls = toc
b_ls = Psilms*a;

fprintf('(|br|, |br-br_ls|^2)=(%4.2e, %4.2e)', norm(b, 2), norm(b-b_ls, 2)^2);
disp('... done');


% get alm from a
alm = cell( maxL+1 , 1);

for ll = 0 : maxL
    
    ind_ls = find( l_grid == ll);
    siz_ls = numel( ind_ls );   
    al = reshape( a(  tablelms(1,:) == ll ), siz_ls, 2*ll+1);
    
    % A_{l,-m} = (-1)^(l+m)* conj( A_{lm} ), thus A_l0 = (-1)^l*conj( A_l0)
    if mod(ll, 2) == 0
        al(:, ll+1) = real( al(:, ll+1) );
    else
        al(:, ll+1) = sqrt(-1)*imag( al(:, ll+1) );
    end
    
    alm{ ll +1 } = al;
end

%% plot Alm(r)
figure(5), clf,
%xx = 0.01: r_mask/100 :r_mask;
%xx = (0:0.1:N)/N;
xx = (0:N-1)/N;
f = zeros( 1, numel( xx));

for ll = 0: 1 : maxL
     
    [X, ~, dems]= xlm( ll, -ll:1:ll, pi/2);
    X = X';
    Y = exp( sqrt(-1)* pi/2* (-ll:1:ll));
    Y = Y.*X;
   for mm = -ll:ll  % nonnegative m, A_{l,-m} = (-1)^(l+m)* conj( A_{lm} )
            
       a_llmm = alm{ ll+1}(:, mm+ll+1);
       
       ind_ls = find( l_grid == ll);
       siz_ls = numel( ind_ls);
       if size( a_llmm, 1) ~= siz_ls
           error('size of a_llmm not equal to siz_ls');
       end
              
       for ss = 1: siz_ls
           root_ls = R_ls( ind_ls( ss )); % s-th root of j_l(r)
           normalcont_ls = sqrt(2)/abs( jl_func(ll+1, root_ls) );%/sqrt(r_mask^3);
           jls = jl_func( ll, root_ls*(xx/r_mask))*normalcont_ls;
           jls = jls*Y(mm+ll+1);
           if ll == 0
               jls(1) = normalcont_ls*Y(1);
           else
               jls(1) = 0;
           end
           f = f + a_llmm(ss)*jls;
       end       
   end
end

%F_a = Psilms(jjorigin:jjorigin+N-1,:)*a;
set(gca, 'FontSize', 16);
hold on;
%plot( xx, real(f),'o-',xx, imag(f), '--g', xx, F_r','+-');
plot( xx, real(f),'o-',xx, imag(f), '--g', xx, b(jjorigin:jjorigin+N-1),'+-');
title(sprintf('Least square interpolation for radial with maxL=%d', ll));
legend('approximation (real part)','approximation (complex part)','original function')
xlabel('radius R/Rmask');
ylabel('function value');

%%
save('cleanrib_gdls_32.mat');
