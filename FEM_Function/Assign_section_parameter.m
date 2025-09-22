%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FE=Assign_section_parameter(FE,Material)
Key=Keywords();
for fi=1:1:numel(FE)
    mat=Material(FE(fi).materialID);
    if strcmp(FE(fi).category,Key.structure)
        FE(fi).detJac=FE(fi).detJac*mat.section.parameter(1);
    end
end
end