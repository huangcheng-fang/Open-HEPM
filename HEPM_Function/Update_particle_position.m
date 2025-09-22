%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FPcell=Update_particle_position(FPcell,FEcell,DOF,U)
CParticles=FPcell.Particles;
CSmoothing_matrix=FPcell.Smoothing_matrix;
CPMid=FPcell.Material_id;
CN_matrix=FEcell.N_matrix;
CElements=FEcell.Elements;
Cdof_num=FEcell.Field_dof_num;
CEMid=FEcell.Material_id;
FEcell_num=numel(CN_matrix);
%==========================================================================
CU=coder.nullcopy(cell(FEcell_num,1));
for fi=1:1:FEcell_num
    numdof=Cdof_num{fi}(1);
    NU=U(DOF(1:numdof,:));
    EU=reshape(NU(:,CElements{fi}),numdof,size(CElements{fi},2),1,[]);
    G_P_U=pagemtimes(EU,'none',CN_matrix{fi},'transpose');
    CU{fi}=reshape(G_P_U,numdof,[]).';
end
%==========================================================================
CU=Collect_cell(CU,CEMid,CPMid,1);
for pci=1:1:numel(CSmoothing_matrix)
    CParticles{pci}=CParticles{pci}+CSmoothing_matrix{pci}*CU{pci};
end
%==============================Output======================================
FPcell.Particles=CParticles;
end