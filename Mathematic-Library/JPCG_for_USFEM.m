%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [x,iter,resF]=JPCG_for_USFEM(K,F,LM,tol,maxiter,x,refF,constraint)
if isempty(LM)
    LM = diag(sparse(1./sum(K.^2,2)));
end
if isempty(x)
    x=zeros(numel(F),1);
end
if isempty(refF)||refF==0
    refF=norm(F);
end

if tol<1e-12
    warning('The tolerance has been modified to 1e-12 for better convergence')
end
if isempty(constraint)
    constraint=zeros(0,2);
end
%==========================================================================
loc=constraint(:,1);
val=constraint(:,2);
PU=zeros(size(K,1),1);
PU(loc,1)=val;
F=F-K*(K'*PU);
if refF==0;refF=norm(F);end
%==========================================================================
Fint=K*(K'*x);
r=F-Fint;
r(loc)=0;
%==========================================================================
resF=norm(r)/refF;iter=0;
if resF<tol;warning('The residual is already less than the tolerance, no iteration required');return;end
%==========================================================================
h=LM*r;
p=h;
%==========================================================================
while resF>tol
    KXP=K*(K'*p);
    alpha=h'*r/(p'*KXP);
    dx=alpha*p;
    x=x+dx;
    beta=h'*r;
    %---------------
    r=r-alpha*KXP;
    r(loc)=0;
    %-----------------
    h=LM*r;
    beta=(h'*r)/beta;
    p=h+beta*p;
    %-----------------
    resF=norm(r)/refF;
    iter=iter+1;
    if iter>maxiter
        error('Exceed maximum number of iterations')
    end
end
x(loc)=val;
%==========================================================================
end