%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [x,iter,resF]=ICCG(K,F,LM,tol,maxiter,x,refF,x_ref,tolscale)
if isempty(LM)
    LM = Get_ichol(K);
end
if isempty(x)
    x=zeros(numel(F),1);
end
if isempty(refF)||refF==0
    refF=norm(F);
end
if isempty(x_ref)
    x_ref=zeros(numel(F),1);
end
if isempty(tolscale)
    tolscale=10;
end
if numel(x)~=size(K,2)||numel(F)~=size(K,2)
    error('Vector and matrix dimensions do not match')
end
if tol<1e-12
    warning('The tolerance has been modified to 1e-12 for better convergence')
end
%==========================================================================
Fint=K*x;
r=F-Fint;
% ref=norm(r);
%==========================================================================
resF=norm(r)/refF,iter=0;
tol=max(min(1e-3,resF/tolscale),tol);
if resF<tol;warning('The residual is already less than the tolerance, no iteration required');return;end
%==========================================================================
UM=LM.';
h=UM\(LM\r);
p=h;
resU=1;
%==========================================================================
while resF>tol||resU>tol
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
    resF=norm(r)/refF;
    resU=norm(dx)/norm(x+x_ref)/10;
    iter=iter+1;
    if iter>maxiter
        error('Exceed maximum number of iterations')
    end
end
%==========================================================================
end

function LM=Get_ichol(K)
try
    LM = ichol(K);
catch
    diagcomp=[1e-3,1e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
    while 1;try LM = ichol(K, struct('diagcomp',diagcomp(1)));break;catch;diagcomp(1)=[];end;end
end
end