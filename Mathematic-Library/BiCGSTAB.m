%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [x,iter,resF]=BiCGSTAB(K,F,LM,tol,maxiter,x)
if isempty(LM)
    [LM,UM] = ilu(K);
end
if isempty(x)
    x=zeros(numel(F),1);
end
if numel(x)~=size(K,2)||numel(F)~=size(K,2)
    error('Vector and matrix dimensions do not match')
end
if tol<1e-12
    warning('The tolerance has been modified to 1e-12 for better convergence')
end
%==========================================================================
refF=norm(F);
Fint=K*x;
r=F-Fint;
rbar=r;
row=rbar'*r;
p=r;iter=0;
%==========================================================================
while resF>tol
    y=UM\(LM\p);
    v=K*y;
    alpha=row/(rbar'*v);
    h=x+alhap*y;
    s=r-alhap*v;
    z=UM\(LM\s);
    t=A*z;
    temp=LM\t;
    w=temp'*(LM\s)/(temp'*temp);
    x=h+w*z;
    r=s-w*t;
    row0=row;
    row=rbar'*r;
    beta=row/row0*(alpha/w);
    p=r+beta*(p-w*v);
    resF=norm(r)/refF;
    iter=iter+1;
    if iter>maxiter
        error('Exceed maximum number of iterations')
    end
end
%==========================================================================
end