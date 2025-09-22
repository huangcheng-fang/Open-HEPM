%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FPcell=Calculate_Kernel_matrix(FPcell,FEcell)
tic;
CDetJac=FEcell.DetJac;
CN_matrix=FEcell.N_matrix;
CNodes=FEcell.Nodes;
CMid=FEcell.Material_id;
FEnum=numel(CMid);

CPMid=FPcell.Material_id;
CParticles=FPcell.Particles;
CInteraction_matrix=FPcell.Interaction_matrix;
CSmoothing_matrix=FPcell.Smoothing_matrix;
CNode_volume=FPcell.Node_volume;
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
%==============================Kernel matrix===============================
for pci=1:1:FPnum
    old_particles=CParticles{pci};
    gauss_points=CGauss_points{pci};
    gauss_volumes=CGauss_volumes{pci};
    r=(gauss_volumes*0.75/pi).^(1/3);
    %----------------------------------------------------------------------
    KDTree=KDTreeSearcher(old_particles);
    [~,Distance] = knnsearch(KDTree,gauss_points,'k',1);
    if ~isempty(Distance)
        new_particles=gauss_points(Distance>=2*r,:);
    else
        new_particles=gauss_points;
    end
    if ~isempty(new_particles)
        warning('new_particles')
        figure(100)
        plot3(new_particles(:,1),new_particles(:,2),new_particles(:,3),'.','MarkerSize',20)
    end
    %----------------------------------------------------------------------
    particles=[old_particles;new_particles];
    KDTree=KDTreeSearcher(particles);
    [PID,Distance] = knnsearch(KDTree,gauss_points,'k',100);
    Distance=r.\Distance;
    %----------------------------------------------------------------------
    Distance_weight=zeros(size(Distance,1),size(Distance,2));
    locD1=Distance<1;
    Distance_weight(locD1)=(4-6*Distance(locD1).^2+3*Distance(locD1).^3);
    locD2=Distance<2&Distance>=1;
    Distance_weight(locD2)=(2-Distance(locD2)).^3;
    locD=Distance<2;
    Distance_weight=sum(Distance_weight,2).\Distance_weight;
    EID=repmat((1:1:size(PID,1))',1,size(PID,2));
    Interaction_matrix=sparse(PID(locD),EID(locD),Distance_weight(locD),size(particles,1),size(PID,1));
    Smoothing_matrix=Interaction_matrix*diag(sparse(gauss_volumes));
    %----------------------------------------------------------------------
    node_volume=full(sum(Smoothing_matrix,2));
    loc=node_volume<=1e-10*mean(node_volume);
    node_volume(loc)=1;
    Smoothing_matrix=diag(sparse(node_volume))\Smoothing_matrix;
    node_volume(loc)=0;
    Pactivation=~loc;
    Smoothing_matrix=diag(sparse(Pactivation))*Smoothing_matrix;
    %----------------------------------------------------------------------
    if abs(1-sum(node_volume,'all')/sum(gauss_volumes,'all'))>1e-8
    warning('gauss_volumes~=node_volume')
    end
    %----------------------------------------------------------------------
    CParticles{pci}=particles;
    CInteraction_matrix{pci}=Interaction_matrix;
    CSmoothing_matrix{pci}=Smoothing_matrix;
    CNode_volume{pci}=node_volume;
end
%===============================Output=====================================
FPcell.Particles=CParticles;
FPcell.Interaction_matrix=CInteraction_matrix;
FPcell.Smoothing_matrix=CSmoothing_matrix;
FPcell.Node_volume=CNode_volume;
%==========================================================================
time=toc;fprintf("Particle Kernel matrix is updated:%fs\n",time);
end