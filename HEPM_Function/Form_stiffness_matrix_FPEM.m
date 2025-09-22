%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function K=Form_stiffness_matrix_FPEM(FEcell,FEUdof,FPcell,Dmatrix,maxdof)
tic
CB_matrix=FEcell.B_matrix;
CEid=FEcell.Element_id;
CEMaterial_id=FEcell.Material_id;
FEcell_num=numel(CEid);
CSmoothing_matrix=FPcell.Smoothing_matrix;
CNode_volume=FPcell.Node_volume;
CPMaterial_id=FPcell.Material_id;
FPcell_num=numel(CPMaterial_id);
for fi=1:FEcell_num
    Bsize=size(CB_matrix{fi});
    KValue=zeros(9,Bsize(2),Bsize(3),Bsize(4));
    dN=[CB_matrix{fi}(1,1:3:end,:,:);CB_matrix{fi}(2,2:3:end,:,:);CB_matrix{fi}(3,3:3:end,:,:);];
    KValue(1:3,1:3:end,:,:)=dN;
    KValue(4:6,2:3:end,:,:)=dN;
    KValue(7:9,3:3:end,:,:)=dN;
    CB_matrix{fi}=KValue;
end
%=========================Calculate_B_matrix===============================
CB_sparse=coder.nullcopy(cell(FPcell_num,1));
for pci=1:1:FPcell_num
    enum=size(CSmoothing_matrix{pci},2);
    B=sparse(enum*9,maxdof);flag=0;
    for fi=1:FEcell_num
        if CEMaterial_id{fi}~=CPMaterial_id{pci};continue;end
        KValue=CB_matrix{fi};
        Ksize=size(KValue,[1,2,3,4]);
        falgend=flag+Ksize(1)*Ksize(3)*Ksize(4);
        KIndexI=reshape(flag+1:falgend,Ksize(1),1,Ksize(3),Ksize(4));
        KIndexI=repmat(KIndexI,1,Ksize(2),1,1);
        KIndexJ=repmat(pagetranspose(FEUdof{fi}),Ksize(1),1,Ksize(3),1);
        flag=falgend;
        B=B+sparse(KIndexI(:),KIndexJ(:),KValue(:),enum*9,maxdof);
    end
    CB_sparse{pci}=B;
end
%======================Calculate_stiffness_matrix==========================
K=sparse(maxdof,maxdof);
for pci=1:1:FPcell_num
    % Form_D_matrix
    NodeDofList2=reshape(1:size(Dmatrix{pci},4)*9,9,[]);
    DV=pagemtimes(Dmatrix{pci},permute(CNode_volume{pci},[4,2,3,1]));
    DI=permute(NodeDofList2,[1,3,2]);
    DI=repmat(DI,1,size(DI,1));
    DJ=pagetranspose(DI);
    D=sparse(DI(:),DJ(:),DV(:));

    % Calculate_K_matrix
    Smoothing_matrix=kron(CSmoothing_matrix{pci},speye(9));
    SB=Smoothing_matrix*CB_sparse{pci};
    K=K+SB'*D*SB;
end
time=toc;fprintf("FPEM stiffness matrix dstrain is calculated:%fs\n",time);
end
