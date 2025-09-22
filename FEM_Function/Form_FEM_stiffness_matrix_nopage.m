%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function K=Form_FEM_stiffness_matrix_nopage(FE,EDep)
tic;
Dof=FE.dof;
DetJac=FE.detJac;
B_matrix=FE.B_matrix;
EIPnum=size(B_matrix,3);
Enum=size(B_matrix,4);
%========================Form_stiffness_matrix========================
dofdim=size(B_matrix,2);
KValue=zeros(dofdim,dofdim,EIPnum,Enum);
KIndexI=zeros(dofdim,dofdim,EIPnum,Enum);
for ei=1:1:Enum
    for ipi=1:1:EIPnum
        KValue(:,:,ei)=KValue(:,:,ei)+B_matrix(:,:,ipi,ei).'*(EDep(:,:,ipi,ei)*DetJac(:,:,ei))*B_matrix(:,:,ipi,ei);
    end
    KIndexI(:,:,ei)=repmat(Dof(:,ei),1,dofdim);
end
KIndexJ=pagetranspose(KIndexI);
%===============================Output=====================================
K.Value=KValue(:);K.IndexI=KIndexI(:);K.IndexJ=KIndexJ(:);
%--------------------------------------------------------------------------
time=toc;fprintf("System stiffness matrix is calculated: %fs\n",time);
end