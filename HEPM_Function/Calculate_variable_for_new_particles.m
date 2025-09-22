%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [FPcell,CPfield]=Calculate_variable_for_new_particles(FPcell,CPfield)
CInteraction_matrix=FPcell.Interaction_matrix;
CSmoothing_matrix=FPcell.Smoothing_matrix;
CNode_volume=FPcell.Node_volume;
CParticles=FPcell.Particles;
% CU=CPfield.U;
CStress=CPfield.Stress;
CStrain=CPfield.Strain;
CPstrain=CPfield.Pstrain;
CParameter=CPfield.Parameter;
global NEfield_G
%==========================================================================
for pci=1:1:numel(CNode_volume)
Interaction_matrix=CInteraction_matrix{pci};
Smoothing_matrix=CSmoothing_matrix{pci};
node_volume=CNode_volume{pci};
particles=CParticles{pci};
%--------------------------------------------------------------------------
Pactivation=node_volume~=0;
Smoothing_matrix=Smoothing_matrix.';
Smoothing_matrix=Smoothing_matrix(:,Pactivation).';
Interaction_matrix=Interaction_matrix.';
Interaction_matrix=Interaction_matrix(:,Pactivation).';
node_volume=node_volume(Pactivation,:);
%--------------------------------------------------------------------------
old_num=size(CStress{pci},4);
particles_old=particles(1:old_num,:);
particles_new=particles(old_num+1:end,:);
Pactivation_old=Pactivation(1:old_num,:);
Pactivation_new=Pactivation(old_num+1:end,:);
particles_new=particles_new(Pactivation_new,:);
particles_old=particles_old(Pactivation_old,:);
particles=[particles_old;particles_new];
%==========================================================================
% U=permute(CU{pci}(:,:,:,Pactivation_old),[1,4,3,2]);
Stress=permute(CStress{pci}(:,:,:,Pactivation_old),[1,4,3,2]);
Strain=permute(CStrain{pci}(:,:,:,Pactivation_old),[1,4,3,2]);
Pstrain=permute(CPstrain{pci}(:,:,:,Pactivation_old),[1,4,3,2]);
Parameter=permute(CParameter{pci}(:,:,:,Pactivation_old),[1,4,3,2]);
%--------------------------------------------------------------------------
% particles_old=NEfield_G.nodes;
KDTree=KDTreeSearcher(particles_old);
[PID1,Distance] = knnsearch(KDTree,particles_new,'k',5);
weight=1./Distance;weight=sum(weight,2).\weight;
PID2=repmat((1:1:size(PID1,1))',1,size(PID1,2));
interp_matrix=sparse(PID1,PID2,weight,size(particles_old,1),size(PID1,1));
%--------------------------------------------------------------------------
% new_U=U*interp_matrix;
new_stress=Stress*interp_matrix;
new_strain=Strain*interp_matrix;
new_pstrain=Pstrain*interp_matrix;
new_parameter=Parameter*interp_matrix;
%==========================================================================
CParticles{pci}=particles;
CInteraction_matrix{pci}=Interaction_matrix;
CSmoothing_matrix{pci}=Smoothing_matrix;
CNode_volume{pci}=node_volume;
% CU{pci}=permute([U,new_U],[1,4,3,2]);
CStress{pci}=permute([Stress,new_stress],[1,4,3,2]);
CStrain{pci}=permute([Strain,new_strain],[1,4,3,2]);
CPstrain{pci}=permute([Pstrain,new_pstrain],[1,4,3,2]);
CParameter{pci}=permute([Parameter,new_parameter],[1,4,3,2]);
end
%==========================================================================
FPcell.Interaction_matrix=CInteraction_matrix;
FPcell.Smoothing_matrix=CSmoothing_matrix;
FPcell.Node_volume=CNode_volume;
FPcell.Particles=CParticles;
% CPfield.U=CU;
CPfield.Stress=CStress;
CPfield.Strain=CStrain;
CPfield.Pstrain=CPstrain;
CPfield.Parameter=CParameter;
end