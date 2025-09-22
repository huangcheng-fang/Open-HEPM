%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [KF,KCP]=Form_FEM_fluid_matrix(FE,maxdof)
tic;
for fi=numel(FE):-1:1
    SDof{fi}=FE(fi).dof;
    SB_matrix{fi}=FE(fi).B_matrix;
    DetJac{fi}=FE(fi).detJac;

    FDof{fi}=SDof{fi}(3:3:end,:)/3;
    FB_matrix{fi}=FE(fi).FB_matrix;

    EDC{fi}=FE(fi).D_conductivity;
    N_matrix{fi}=FE(fi).N_matrix;
end
%===========================Form_stiffness_matrix==========================
KFs=Create_parallel_container(numel(DetJac));
for fi=1:1:numel(DetJac)
    if isempty(EDC{fi});continue;end
    EDepJac=pagemtimes(EDC{fi},DetJac{fi});
    KValue=pagemtimes(FB_matrix{fi},'transpose',EDepJac,'none');
    KValue=pagemtimes(KValue,FB_matrix{fi});
    KValue=permute(sum(KValue,3),[1,2,4,3]);
    %--------------------------------------------------------------------------
    KIndexI=permute(FDof{fi},[1,3,2]);
    KIndexI=repmat(KIndexI,1,size(KIndexI,1),1);
    KIndexJ=pagetranspose(KIndexI);
    KFs(fi,1).c1=KIndexI(:);
    KFs(fi,1).c2=KIndexJ(:);
    KFs(fi,1).c3=KValue(:);
end
[KFI,KFJ,KFV]=Merge_container(KFs,1);
%===========================Form_coupling_matrix===========================
KCPs=Create_parallel_container(numel(DetJac));
for fi=1:1:numel(DetJac)
    if isempty(EDC{fi});continue;end
    CValue=pagemtimes(DetJac{fi},N_matrix{fi});
    CValue=pagemtimes(sum(SB_matrix{fi}(1:3,:,:,:),1),'transpose',CValue,'none');
    CValue=permute(sum(CValue,3),[1,2,4,3]);
    %--------------------------------------------------------------------------
    CIndexI=permute(SDof{fi},[1,3,2]);
    CIndexJ=permute(FDof{fi},[3,1,2]);
    CIndexI=repmat(CIndexI,1,size(CIndexJ,2),1);
    CIndexJ=repmat(CIndexJ,size(CIndexI,1),1,1);
   
    KCPs(fi,1).c1=CIndexI(:);
    KCPs(fi,1).c2=CIndexJ(:);
    KCPs(fi,1).c3=CValue(:);
end
[KCPI,KCPJ,KCPV]=Merge_container(KCPs,1);
%===============================Output=====================================
KF=sparse(KFI,KFJ,KFV,maxdof/3,maxdof/3);
KCP=sparse(KCPI,KCPJ,KCPV,maxdof,maxdof/3);
%--------------------------------------------------------------------------
time=toc;fprintf("FEM Fluid stiffness matrix is calculated: %fs\n",time);
end