%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function M=Form_mass_vector(FE,Material,maxdof)
tic;
%==============================Input=======================================
for fi=numel(FE):-1:1
    DetJac{fi}=FE(fi).detJac;
    Dof{fi}=FE(fi).dof;
    N_matrix{fi}=FE(fi).N_matrix;
    density{fi}=Material(FE(fi).materialID).elasticity.parameter(3);
end
%========================Form_mass_vector========================
Mcontainer=Create_parallel_container(numel(DetJac));
for fi=1:1:numel(DetJac)
    DetJacdensity=pagemtimes(DetJac{fi},density{fi});
    MValue=pagemtimes(N_matrix{fi},'transpose',DetJacdensity,'none');
    MValue=repmat(permute(sum(MValue,3),[1,2,4,3]),1,3,1);
    %----------------------------------------------------------------------
    MIndexI=permute(Dof{fi},[1,3,2]);
    MIndexI=reshape(MIndexI,3,[],size(MIndexI,3));
    MIndexI=pagetranspose(MIndexI);
    Mcontainer(fi,1).c1=MIndexI(:);
    Mcontainer(fi,1).c2=MValue(:);
end
[I,V]=Merge_container(Mcontainer,1);
%===============================Output=====================================
M=accumarray(I,V,[maxdof,1]);
%--------------------------------------------------------------------------
time=toc;fprintf("Mass vector is calculated: %fs\n",time);
end