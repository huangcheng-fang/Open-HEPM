%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [x,iter,resF]=PCG_IC(K,F,LM,tol,maxiter,x0)
if nargin<6
    x0=zeros(numel(F),1);
end
if numel(x0)~=size(K,2)||numel(F)~=size(K,2)
    error('Vector and matrix dimensions do not match')
end
%==========================================================================
ref=norm(F);
UM=LM.';
x=x0;
%==========================================================================
r=F-K*x;
h=UM\(LM\r);
p=h;
resF=norm(r)/ref;iter=0;
%==========================================================================
while resF>tol
    KXP=K*p;
    alpha=h'*r/(p'*KXP);
    dx=alpha*p;
    x=x+dx;
    beta=h'*r;
    %---------------
    r=r-alpha*KXP;
    %-----------------
    h=UM\(LM\r);
    beta=(h'*r)/beta;
    p=h+beta*p;
    %-----------------
    resF=norm(r)/ref;
    
    iter=iter+1;
    if iter>maxiter
        error('Exceed maximum number of iterations')
    end
end
%==========================================================================
end