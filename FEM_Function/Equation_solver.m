%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [X,iter]=Equation_solver(K,F,tol,type)
tic
%--------------------------------------------------------------------------
if strcmp(type,'default')
    if ~ismember(0,diag(K))&&full(sum(abs(K.'-K),'all')/sum(abs(K.'),'all'))<1e-8&&isempty(find(diag(K)<0, 1))
        type='IC-CG';
    elseif min(diag(K))>0
        type='J-CG(nonsymmetry)';
    elseif min(diag(K))==0
        type='ILU-BCG';
    else
        type='Direct';
    end
end
%--------------------------------------------------------------------------
maxiter=2000;
D=diag(K);D(D==0)=1;D=sparse(1./sqrt(D));D=diag(D);K=D*K*D;F=D*F;
switch type
    case 'IC-CG'
        try
            LM = ichol(K);UM=LM.';
        catch
            diagcomp=[1e-3,1e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
            while 1;try LM = ichol(K, struct('diagcomp',diagcomp(1)));break;catch;diagcomp(1)=[];end;end
            UM=LM.';type=['IC-CG(diagcomp=',num2str(diagcomp(1)),')'];
        end
        [X,flag,relres,iter]= pcg(K,F,tol,maxiter,LM,UM);
    case 'J-CG(nonsymmetry)'
        [X,flag,relres,iter] = pcg(K,F,tol,maxiter);
    case 'ILU-BCG'
        [LM,UM] = ilu(K);
        [X,flag,relres,iter]= bicgstab(K,F,tol,maxiter,LM,UM);
    case 'Direct'
        X=K\F;relres=0;flag=0;type='Direct-Method';iter=1;
    otherwise
        error('Unkown equation solver type')
end
X=D*X;
%==========================================================================
if flag||relres>tol;warning(['Iterative non-convergence:',num2str(relres)]);end
time=toc;disp(['Iteration number of ',type,' is:',num2str(iter),'(time=',num2str(time),'s)'])
end