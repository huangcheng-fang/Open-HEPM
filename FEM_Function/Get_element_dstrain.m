%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [CEdstrain,CEspin]=Get_element_dstrain(stepU,FEcell,FEUdof,GNL)
CType=FEcell.Type;
CNodes=FEcell.Nodes;
CB_matrix=FEcell.B_matrix;
CDetJac=FEcell.DetJac;
NumFEcell=numel(CType);
%===========================initialization=================================
CEdstrain=cell(NumFEcell,1);
CEspin=cell(NumFEcell,1);
%============================small_strain==================================
if ~GNL
    for fi=1:1:NumFEcell
        udof=FEUdof{fi};
        EdU=stepU(udof);
        dstrain=pagemtimes(CB_matrix{fi},EdU);
        CEdstrain{fi}=dstrain;
    end
end
%============================large_strain==================================
if GNL
    for fi=1:1:NumFEcell
        udof=FEUdof{fi};
        CNodes_stepU=pagetranspose(reshape(stepU(udof),3,size(CNodes{fi},1),1,[]));
        switch CType{fi}
            case 334
                [CEdstrain{fi},CEspin{fi}]=Tetrahedron_large_strain(CNodes{fi},CNodes_stepU);
            case 338
                IntegralOrder=size(CB_matrix{fi},3)^(1/3);
                [CEdstrain{fi},CEspin{fi}]=Hexahedron_large_strain(CNodes{fi},CNodes_stepU,IntegralOrder);
            case 3381
                IntegralOrder=size(CB_matrix{fi},3)^(1/3);
                [CEdstrain{fi},CEspin{fi}]=Hexahedron_large_strain(CNodes{fi},CNodes_stepU,IntegralOrder);
                Ev=mean(CEdstrain{fi}(1:3,:,:,:),1);
                Ev_mean=sum(pagemtimes(Ev,CDetJac{fi}),3)./sum(CDetJac{fi},3);
                CEdstrain{fi}(1:3,:,:,:)=CEdstrain{fi}(1:3,:,:,:)-repmat(Ev,3,1,1)+repmat(Ev_mean,3,1,size(Ev,3));
            otherwise
                error('Unknown element type in large_strain')
        end
    end
end
end