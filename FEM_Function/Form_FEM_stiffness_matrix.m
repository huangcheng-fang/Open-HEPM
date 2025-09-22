%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function K=Form_FEM_stiffness_matrix(FEcell,FEUdof,Dmatrix,maxdof)
tic;
CDetJac=FEcell.DetJac;
CB_matrix=FEcell.B_matrix;
NumFEcell=numel(CDetJac);
for fi=1:NumFEcell
    Bsize=size(CB_matrix{fi});
    KValue=zeros(9,Bsize(2),Bsize(3),Bsize(4));
    dN=[CB_matrix{fi}(1,1:3:end,:,:);CB_matrix{fi}(2,2:3:end,:,:);CB_matrix{fi}(3,3:3:end,:,:);];
    KValue(1:3,1:3:end,:,:)=dN;
    KValue(4:6,2:3:end,:,:)=dN;
    KValue(7:9,3:3:end,:,:)=dN;
    CB_matrix{fi}=KValue;
end
%=========================Form_stiffness_matrix============================
K=sparse(maxdof,maxdof);
for fi=1:NumFEcell
    EDepJac=pagemtimes(Dmatrix{fi},CDetJac{fi});
    KValue=pagemtimes(CB_matrix{fi},'transpose',EDepJac,'none');
    KValue=pagemtimes(KValue,CB_matrix{fi});
    KValue=permute(sum(KValue,3),[1,2,4,3]);
    %----------------------------------------------------------------------
    KIndexI=FEUdof{fi};
    KIndexI=repmat(KIndexI,1,size(KIndexI,1),1);
    KIndexJ=pagetranspose(KIndexI);
    %----------------------------------------------------------------------
    K=K+sparse(KIndexI(:),KIndexJ(:),KValue(:),maxdof,maxdof);
end
%==========================================================================
time=toc;fprintf("System stiffness matrix is calculated: %fs\n",time);
end