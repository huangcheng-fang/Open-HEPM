%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function C=Form_damping_vector(FE,Material,M,alpha)
tic;
%==============================Input=======================================
for fi=numel(FE):-1:1
    Dof{fi}=FE(fi).dof;
    density=Material(FE(fi).materialID).elasticity.parameter(3);
    E=Material(FE(fi).materialID).elasticity.parameter(1);
    damping_coeff{fi}=alpha*sqrt(E)*density^(-1/6);
end
%========================Form_damping_vector========================
C=zeros(size(M,1),1);
for fi=1:1:numel(FE)
    MIndexI=Dof{fi}(:);
    C(MIndexI)=damping_coeff{fi};
end
%===============================Output=====================================
C=C.*(M).^(2/3);
%--------------------------------------------------------------------------
time=toc;fprintf("Damping Vector is calculated: %fs\n",time);
end