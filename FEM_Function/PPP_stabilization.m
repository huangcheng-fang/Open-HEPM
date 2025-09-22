%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function KF=PPP_stabilization(FE,Material,maxdof)
tic;
%==============================Input=======================================
for fi=numel(FE):-1:1
    type{fi}=FE(fi).type;
    DetJac{fi}=FE(fi).detJac;
    FDof{fi}=FE(fi).dof(3:3:end,:)/3;
    EDC{fi}=FE(fi).D_conductivity;
    N_matrix{fi}=FE(fi).N_matrix;
    pm=Material(FE(fi).materialID).elasticity.parameter(1:2);
    PPPCoef=Material(FE(fi).materialID).fluid.stabilization;
    G=(pm(1)/2/(1+pm(2)));
    PPPfactor{fi}=PPPCoef./G/2;
    if G==0;PPPfactor{fi}=0;end
end
%========================Form_PPP_matrix========================
KFs=Create_parallel_container(numel(DetJac));
for fi=1:1:numel(DetJac)
    if isempty(EDC{fi});continue;end
    H=PPP_matrix(type{fi},N_matrix{fi});
    KValue=pagemtimes(DetJac{fi},PPPfactor{fi});
    KValue=pagemtimes(KValue,H);
    KValue=permute(sum(KValue,3),[1,2,4,3]);
    %----------------------------------------------------------------------
    KIndexI=permute(FDof{fi},[1,3,2]);
    KIndexI=repmat(KIndexI,1,size(KIndexI,1),1);
    KIndexJ=pagetranspose(KIndexI);
    KFs(fi,1).c1=KIndexI(:);
    KFs(fi,1).c2=KIndexJ(:);
    KFs(fi,1).c3=KValue(:);
end
[KFI,KFJ,KFV]=Merge_container(KFs,1);
%===============================Output=====================================
KF=sparse(KFI,KFJ,KFV,maxdof/3,maxdof/3);
%--------------------------------------------------------------------------
time=toc;fprintf("FEM Fluid PPP matrix is calculated: %fs\n",time);
end

function H=PPP_matrix(type,N)
switch type
    case 334
        a=3/80;b=-1/80;
        H=[a,b,b,b;b,a,b,b;b,b,a,b;b,b,b,a];
    case 323
        a=1/18;b=-1/36;
        H=[a,b,b;b,a,b;b,b,a];
    case {338,3381}
        N=N-0.125;
        H=pagemtimes(N,'transpose',N,'none');
    case 312
        N=N-0.5;
        H=pagemtimes(N,'transpose',N,'none');
    case 324
        N=N-0.25;
        H=pagemtimes(N,'transpose',N,'none');
    otherwise
        H=[];
        error('Unknown element type')
end
end