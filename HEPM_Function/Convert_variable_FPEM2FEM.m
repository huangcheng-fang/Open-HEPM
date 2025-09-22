%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Efieldcell=Convert_variable_FPEM2FEM(CPfield,FPcell,CEDetJac,CEMaterial_id,CElement_id)
FEcell_num=numel(CEMaterial_id);
CStress=coder.nullcopy(cell(FEcell_num,1));
CStrain=coder.nullcopy(cell(FEcell_num,1));
CPstrain=coder.nullcopy(cell(FEcell_num,1));
CParameter=coder.nullcopy(cell(FEcell_num,1));
CSmoothing_matrix=FPcell.Smoothing_matrix;
CNode_volume=FPcell.Node_volume;
CPMaterial_id=FPcell.Material_id;
FPcell_num=numel(CPMaterial_id);
%==========================================================================
for pci=1:1:FPcell_num
    node_volume=CNode_volume{pci};
    Smoothing_matrix=CSmoothing_matrix{pci};
    convert_matrix=diag(sparse(node_volume))*Smoothing_matrix;
    Estress=permute(CPfield.Stress{pci},[1,4,3,2])*convert_matrix;Estress=permute(Estress,[1,4,3,2]);
    Estrain=permute(CPfield.Strain{pci},[1,4,3,2])*convert_matrix;Estrain=permute(Estrain,[1,4,3,2]);
    Epstrain=permute(CPfield.Pstrain{pci},[1,4,3,2])*convert_matrix;Epstrain=permute(Epstrain,[1,4,3,2]);
    Eparameter=permute(CPfield.Parameter{pci},[1,4,3,2])*convert_matrix;Eparameter=permute(Eparameter,[1,4,3,2]);
    %==========================================================================
    flag=0;
    for fci=1:1:FEcell_num
        if CEMaterial_id{fci}==CPMaterial_id{pci}
            eid=flag+1:flag+size(CEDetJac{fci},4);
            flag=flag+size(CEDetJac{fci},4);
            CStress{fci}=pagemldivide(CEDetJac{fci},Estress(:,:,:,eid));
            CStrain{fci}=pagemldivide(CEDetJac{fci},Estrain(:,:,:,eid));
            CPstrain{fci}=pagemldivide(CEDetJac{fci},Epstrain(:,:,:,eid));
            CParameter{fci}=pagemldivide(CEDetJac{fci},Eparameter(:,:,:,eid));
        end
    end
end
%=============================Output=======================================
Efieldcell.Stress=CStress;
Efieldcell.Strain=CStrain;
Efieldcell.Pstrain=CPstrain;
Efieldcell.Parameter=CParameter;
Efieldcell.Element_id=CElement_id;
end