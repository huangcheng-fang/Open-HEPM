%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [FPcell,CPfield]=Creat_FPcell(PMID,MID,particles,Pfield)
coder.varsize('Particles')
Particles=coder.nullcopy(cell(1,1));
coder.varsize('Particle_id')
Particle_id=coder.nullcopy(cell(1,1));
coder.varsize('Material_id')
Material_id=coder.nullcopy(cell(1,1));
PMID_list=unique(MID);
%==========================================================================
for mi=1:1:numel(PMID_list)
    pid=find(PMID==PMID_list(mi));
    Particles{mi,1}=particles(pid,:);
    Particle_id{mi,1}=pid;
    Material_id{mi,1}=PMID_list(mi);
end
%==========================================================================
blank=coder.nullcopy(cell(numel(Material_id),1));
FPcell.Particles=Particles;
FPcell.Interaction_matrix=blank;
FPcell.Smoothing_matrix=blank;
FPcell.Node_volume=blank;
FPcell.Material_id=Material_id;
%==========================================================================
% CPfield.U=blank;
CPfield.Stress=blank;
CPfield.Strain=blank;
CPfield.Pstrain=blank;
CPfield.Parameter=blank;
for fi=1:1:numel(Particle_id)
    pid=Particle_id{fi};
    % CPfield.U{fi}=Pfield.U(:,:,:,pid);
    CPfield.Stress{fi}=Pfield.Stress(:,:,:,pid);
    CPfield.Strain{fi}=Pfield.Strain(:,:,:,pid);
    CPfield.Pstrain{fi}=Pfield.Pstrain(:,:,:,pid);
    CPfield.Parameter{fi}=Pfield.Parameter(:,:,:,pid);
end
end