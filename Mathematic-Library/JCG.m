%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [x,iter,res]=JCG(K,F,constraint,tol,maxiter,x0)
if nargin<6
    x0=zeros(numel(F),1);
end
if numel(x0)~=size(K,2)||numel(F)~=size(K,2)
    error('Vector and matrix dimensions do not match')
end
if tol<1e-12
    warning('The tolerance has been modified to 1e-12 for better convergence')
end
%==========================================================================
Fint=K*x0;
Fint(constraint(:,1))=0;
%==========================================================================
PU=zeros(size(K,1),1);
PU(constraint(:,1),1)=constraint(:,2);
F=F-K*PU;ref=norm(F);
F(constraint(:,1))=constraint(:,2);
Aeye=sparse(constraint(:,1),constraint(:,1),1,numel(F),numel(F));
Ceye=speye(numel(F))-Aeye;
K=Ceye*K*Ceye;K=K+Aeye;
%==========================================================================
F=F-Fint;
%==========================================================================
LM = diag(K);
x=zeros(numel(F),1);
Fint=zeros(numel(F),1);
%==========================================================================
r=F-Fint;
h=LM.\r;
p=h;
res=norm(r)/ref;iter=0;
tol=max(min(tol,res/10),1e-12);
%==========================================================================
while res>tol
    KXP=K*p;
    alpha=h'*r/(p'*KXP);
    x=x+alpha*p;
    beta=h'*r;
    %---------------
    Fint=Fint+alpha*KXP;
    r=F-Fint;
    %-----------------
    h=LM.\r;
    beta=(h'*r)/beta;
    p=h+beta*p;
    %-----------------
    res=norm(r)/ref;
    iter=iter+1;
    if iter>maxiter
        error('Exceed maximum number of iterations')
    end
end
%==========================================================================
x=x+x0;
end

function LM=Get_ichol(K)
try
    LM = ichol(K);
catch
    diagcomp=[1e-3,1e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
    while 1;try LM = ichol(K, struct('diagcomp',diagcomp(1)));break;catch;diagcomp(1)=[];end;end
end
end