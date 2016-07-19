function [X, err_rhs, iter, resvec, rhsn]=solve_cg_shrink_deb(A,C, C_tr, w, noise_variance, nim, freq, regu)
% solve problem  \argmin_{X} \sum_i w(i) \|A{i} *X*A{i}'-C{i}\|_F^2 such that X is real
D=size(A{1});
E=zeros(D(2),D(2));
E_pop=zeros(D(2),D(2));
C_tot=zeros(size(C{1}));
rhs_tr=zeros(D(2),D(2));
err_rhs=1;
iter=0;
resvec=0;
for i=1:length(A)
    B{i}=A{i}'*A{i};
    E=E+w(i)*A{i}'*C{i}*A{i};
    C_tot=C_tot+w(i)*C{i};
    E_pop=E_pop + w(i)*A{i}'*A{i};
    rhs_tr=rhs_tr+w(i)*B{i}*C_tr*B{i};
end

[UU,SS]=svd(E_pop);
sqrt_Epop=UU*(SS^0.5)*UU';
inv_sqrt_Epop=inv(sqrt_Epop); %% matrix S

E=inv_sqrt_Epop*E*inv_sqrt_Epop;
rhs_tr=inv_sqrt_Epop*rhs_tr*inv_sqrt_Epop;
[UU,EE]=eig(E);
%[p,fre]= check_MP(diag(EE), freq, noise_variance, nim, 20);

if (freq==0)
    gamma=size(E,1)/nim;
else
    gamma=size(E,1)/(2*nim);
end

%
[EE,ids]=sort(diag(EE),'descend');
UU=UU(:,ids);

if noise_variance~=0
[k_hat, nnv]=KN_rankEst(EE,nim,1,0.1,size(EE,1));
top_kn=zeros(size(EE));
top_kn(1:k_hat)=ones(size(top_kn(1:k_hat)));
EE=EE.*top_kn;
end

if noise_variance~=0
    E= UU*diag(op_shrink(diag(EE), noise_variance, gamma))*UU';
end

%[ev,evals]=eig(C_tot/nim);
%[p,fre]= check_MP(diag(evals), freq, noise_variance, nim, 20);

err_rhs=norm((E-rhs_tr),'fro');
rhsn=norm(rhs_tr,'fro');
if norm(E,'fro')<1e-15
    sprintf('Skipping CG tol 1e-15 for k=%d',freq)
    X=zeros(size(E));
    return;
end

if freq==0
    reg_flag=1;
    %regu=choose_regu_k0_cgshrink(A,inv_sqrt_Epop,w);
else
    reg_flag=0;
end

[X, resvec, iter]=solve_equation(B,E,w,inv_sqrt_Epop,reg_flag, regu);% it is equavalent to solve \sum_i w(i) A{i}'A{i}XA{i}'A{i}=\sum_i w(i) A{i}'CA{i} with X being real
% Scaling both LHS and RHS by nim
end

function [X, resvec, iter]=solve_equation(A, C, w, inv_sqrt_Epop_, reg_flag, regu)
% solve equation \sum_i A{i} *X*A{i}'=C
[D1,D2]=size(A{1});% D2 dimension of C, D1 dimension of X
AA=zeros(D1^2,D2^2);

   
	function [lhs]= applyop(x)
        lhs=0;
        x=reshape(x,D1,D1);
        for i=1:length(A)
            lhs=lhs+w(i)*A{i}*inv_sqrt_Epop_*x*inv_sqrt_Epop_*A{i};% To make LHS PSD so that CG can be used
        end
        
        lhs=reshape(lhs,D1^2,1);
        if reg_flag==1
            lhs=lhs+regu*reshape(x,D1^2,1);
        end
        
        lhs=reshape(lhs,D1,D1);
        lhs=inv_sqrt_Epop_*lhs*inv_sqrt_Epop_;
        lhs=reshape(lhs,D1^2,1);
        
	end
tol = 1e-10;
maxit = 50000;

[X,flag2,relres, iter, resvec] = pcg(@applyop,reshape(real(C),D1^2,1),tol,maxit);
%flag2, iter
X=reshape(X,D2,D2);
X=inv_sqrt_Epop_*X*inv_sqrt_Epop_; % Sigma from Sigma_S
end
