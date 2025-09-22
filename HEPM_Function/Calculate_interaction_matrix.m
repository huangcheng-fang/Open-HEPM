%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [FPcell,CPfield]=Calculate_interaction_matrix(FPcell,CPfield,FEcell)
tic;
CDetJac=FEcell.DetJac;
CN_matrix=FEcell.N_matrix;
CNodes=FEcell.Nodes;
CMid=FEcell.Material_id;
FEnum=numel(CMid);

CStress=CPfield.Stress;
CStrain=CPfield.Strain;
CPstrain=CPfield.Pstrain;
CParameter=CPfield.Parameter;

CPMid=FPcell.Material_id;
CParticles=FPcell.Particles;
CInteraction_matrix=FPcell.Interaction_matrix;
FPnum=numel(CPMid);
%=============================gauss point==================================
CGauss_points=coder.nullcopy(cell(FEnum,1));
CGauss_volumes=coder.nullcopy(cell(FEnum,1));
for fci=1:1:FEnum
    ip=reshape(pagemtimes(CN_matrix{fci},CNodes{fci}),3,[]);
    CGauss_points{fci}=ip.';
    CGauss_volumes{fci}=CDetJac{fci}(:);
end
CGauss_points=Collect_cell(CGauss_points,CMid,CPMid,1);
CGauss_volumes=Collect_cell(CGauss_volumes,CMid,CPMid,1);
%=============================add particle==================================
for pci=1:1:FPnum
    particles_old=CParticles{pci};
    gauss_points=CGauss_points{pci};
    gauss_volumes=CGauss_volumes{pci};
    r=(gauss_volumes*0.75/pi).^(1/3);
    %----------------------------------------------------------------------
    KDTree=KDTreeSearcher(particles_old);
    [~,Distance] = knnsearch(KDTree,gauss_points,'k',1);
    if ~isempty(Distance)
        particles_new=gauss_points(Distance>=1.8*r,:);
    else
        particles_new=gauss_points;
    end
    if ~isempty(particles_new)
1
    end
    %----------------------------------------------------------------------
    Stress=permute(CStress{pci},[1,4,3,2]);
    Strain=permute(CStrain{pci},[1,4,3,2]);
    Pstrain=permute(CPstrain{pci},[1,4,3,2]);
    Parameter=permute(CParameter{pci},[1,4,3,2]);
    %----------------------------------------------------------------------
    KDTree=KDTreeSearcher(particles_old);
    [PID1,Distance] = knnsearch(KDTree,particles_new,'k',5);
    weight=1./Distance;weight=sum(weight,2).\weight;
    PID2=repmat((1:1:size(PID1,1))',1,size(PID1,2));
    interp_matrix=sparse(PID1,PID2,weight,size(particles_old,1),size(PID1,1));
    %----------------------------------------------------------------------
    new_stress=Stress*interp_matrix;
    new_strain=Strain*interp_matrix;
    new_pstrain=Pstrain*interp_matrix;
    new_parameter=Parameter*interp_matrix;
    %----------------------------------------------------------------------
    CParticles{pci}=[particles_old;particles_new];
    CStress{pci}=permute([Stress,new_stress],[1,4,3,2]);
    CStrain{pci}=permute([Strain,new_strain],[1,4,3,2]);
    CPstrain{pci}=permute([Pstrain,new_pstrain],[1,4,3,2]);
    CParameter{pci}=permute([Parameter,new_parameter],[1,4,3,2]);
end
%===========================Interaction_matrix=============================
for pci=1:1:FPnum
    particles=CParticles{pci};
    gauss_points=CGauss_points{pci};
    gauss_volumes=CGauss_volumes{pci};
    r=(gauss_volumes*0.75/pi).^(1/3);
    %----------------------------------------------------------------------
    KDTree=KDTreeSearcher(particles);
    [PID,Distance] = knnsearch(KDTree,gauss_points,'k',10);
    Distance=r.\Distance;
    %----------------------------------------------------------------------
    strength=Interaction_strength(Distance);
    %----------------------------------------------------------------------
    loc=strength~=0;
    EID=repmat((1:1:size(PID,1))',1,size(PID,2));
    Interaction_matrix=sparse(PID(loc),EID(loc),strength(loc),size(particles,1),size(PID,1));
    %----------------------------------------------------------------------
    CParticles{pci}=particles;
    CInteraction_matrix{pci}=Interaction_matrix;
end
%===========================remove particle================================
for pci=1:1:FPnum
Interaction_matrix=CInteraction_matrix{pci};
Pactivation=sum(Interaction_matrix,2)>=8e-3;
%--------------------------------------------------------------------------
Interaction_matrix=Interaction_matrix.';
CParticles{pci}=CParticles{pci}(Pactivation,:);
%--------------------------------------------------------------------------
CInteraction_matrix{pci}=Interaction_matrix(:,Pactivation).';
CStress{pci}=CStress{pci}(:,:,:,Pactivation);
CStrain{pci}=CStrain{pci}(:,:,:,Pactivation);
CPstrain{pci}=CPstrain{pci}(:,:,:,Pactivation);
CParameter{pci}=CParameter{pci}(:,:,:,Pactivation);
end
%===============================Output=====================================
FPcell.Particles=CParticles;
FPcell.Interaction_matrix=CInteraction_matrix;
CPfield.Stress=CStress;
CPfield.Strain=CStrain;
CPfield.Pstrain=CPstrain;
CPfield.Parameter=CParameter;
%==========================================================================
time=toc;fprintf("Particle-point interaction matrix is calculated:%fs\n",time);
end


function strength=Interaction_strength(Distance)
strength=zeros(size(Distance,1),size(Distance,2));
locD1=Distance<1;
strength(locD1)=(4-6*Distance(locD1).^2+3*Distance(locD1).^3);
locD2=Distance<2&Distance>=1;
strength(locD2)=(2-Distance(locD2)).^3;
end