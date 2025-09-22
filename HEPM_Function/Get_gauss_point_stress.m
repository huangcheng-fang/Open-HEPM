%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function CEStress=Get_gauss_point_stress(CStress,FPcell,CEDetJac,CEMaterial_id)
FEcell_num=numel(CEMaterial_id);
CEStress=coder.nullcopy(cell(FEcell_num,1));
CSmoothing_matrix=FPcell.Smoothing_matrix;
CNode_volume=FPcell.Node_volume;
CPMaterial_id=FPcell.Material_id;
FPcell_num=numel(CPMaterial_id);
%==========================================================================
for pci=1:1:FPcell_num
    node_volume=CNode_volume{pci};
    Smoothing_matrix=CSmoothing_matrix{pci};
    convert_matrix=diag(sparse(node_volume))*Smoothing_matrix;
    Estress=permute(CStress{pci},[1,4,3,2])*convert_matrix;
    Estress=permute(Estress,[1,4,3,2]);
    flag=0;
    for fci=1:1:FEcell_num
        if CEMaterial_id{fci}==CPMaterial_id{pci}
            eid=flag+1:flag+size(CEDetJac{fci},4);
            flag=flag+size(CEDetJac{fci},4);
            CEStress{fci}=pagemldivide(CEDetJac{fci},Estress(:,:,:,eid));
        end
    end
end
end