%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Material,Contact_condition]=Permeability_X_Time(Material,Contact_condition,dT)
for mi=1:1:numel(Material)
    Material(mi).fluid.conductivity=Material(mi).fluid.conductivity*dT;
end
for ci=1:1:numel(Contact_condition)
    Contact_condition(ci).parameter.permeability=Contact_condition(ci).parameter.permeability*dT;
end
end