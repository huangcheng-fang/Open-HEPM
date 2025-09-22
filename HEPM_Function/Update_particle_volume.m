%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FPcell=Update_particle_volume(FPcell,FEcell)
tic;
CDetJac=FEcell.DetJac;
CMid=FEcell.Material_id;
FEnum=numel(CMid);

CPMid=FPcell.Material_id;
CInteraction_matrix=FPcell.Interaction_matrix;
CNode_volume=FPcell.Node_volume;
FPnum=numel(CPMid);
%=============================gauss point==================================
CGauss_volumes=coder.nullcopy(cell(FEnum,1));
for fci=1:1:FEnum
    CGauss_volumes{fci}=CDetJac{fci}(:);
end
CGauss_volumes=Collect_cell(CGauss_volumes,CMid,CPMid,1);
%==============================Kernel matrix===============================
for pci=1:1:FPnum
    Interaction_matrix=CInteraction_matrix{pci};
    gauss_volumes=CGauss_volumes{pci};
    %----------------------------------------------------------------------
    Smoothing_matrix=Interaction_matrix*diag(sparse(gauss_volumes)./sum(Interaction_matrix,1).');
    %----------------------------------------------------------------------
    node_volume=full(sum(Smoothing_matrix,2));
    Smoothing_matrix=diag(sparse(node_volume))\Smoothing_matrix;
    %----------------------------------------------------------------------
    if abs(1-sum(node_volume,'all')/sum(gauss_volumes,'all'))>1e-8
    warning('gauss_volumes~=node_volume')
    end
    %----------------------------------------------------------------------
    CInteraction_matrix{pci}=Interaction_matrix;
    CSmoothing_matrix{pci}=Smoothing_matrix;
    CNode_volume{pci}=node_volume;
end
%===============================Output=====================================
FPcell.Smoothing_matrix=CSmoothing_matrix;
FPcell.Node_volume=CNode_volume;
%==========================================================================
time=toc;fprintf("Particle volume is updated:%fs\n",time);
end