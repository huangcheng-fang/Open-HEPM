%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function G_vector=Form_FEM_fluid_gravity_vector(FE,maxdof)
tic;
for fi=numel(FE):-1:1
    SDof{fi}=FE(fi).dof;
    DetJac{fi}=FE(fi).detJac;

    FDof{fi}=SDof{fi}(3:3:end,:)/3;
    FB_matrix{fi}=FE(fi).FB_matrix;
    R_matrix{fi}=FE(fi).R_matrix;
    
    EDC{fi}=FE(fi).D_conductivity;
    weight{fi}=FE(fi).fluid_weight(1,:,:,:);
    direct{fi}=FE(fi).fluid_weight(2,1,1,1);
end
%===========================Form_fluid_gravity_vector==========================
G_vectors=Create_parallel_container(numel(DetJac));
for fi=1:1:numel(DetJac)
    if isempty(EDC{fi});continue;end
    if ~isempty(R_matrix{fi})
        Rweight=pagemtimes(R_matrix{fi}(1:size(FB_matrix{fi},1),direct{fi},:,:),weight{fi});
    else
        Rweight=zeros(3,1);
        Rweight(direct{fi})=weight{fi};
    end
    EDepJac=pagemtimes(EDC{fi},Rweight);
    EDepJac=pagemtimes(EDepJac,DetJac{fi});
    GValue=pagemtimes(FB_matrix{fi},'transpose',EDepJac,'none');
    GValue=permute(sum(GValue,3),[1,2,4,3]);
    %--------------------------------------------------------------------------
    GIndexI=permute(FDof{fi},[1,3,2]);
    G_vectors(fi,1).c1=GIndexI(:);
    G_vectors(fi,1).c2=GValue(:);
end
[GI,GV]=Merge_container(G_vectors,1);
%===============================Output=====================================
if isempty(GI)
    G_vector=zeros(maxdof/3,1);
else
    G_vector=accumarray(GI,GV,[maxdof/3,1]);
end
%--------------------------------------------------------------------------
time=toc;fprintf("FEM Fluid gravity vector is calculated: %fs\n",time);
end