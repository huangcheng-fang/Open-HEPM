%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Model,Result]=Mesh_adjustment_solver(Model,Result)
NResult=Mapping_variable_E2N(Result,Model.Mesh);
[Nset,Sset,Eset]=Save_set_information(Model.Set,Model.Mesh);
%--------------------------------------------------------------------------
Model.Mesh=Remesh_domain(Model.Mesh,[]);
%--------------------------------------------------------------------------
Result=Mapping_variable_N2E(NResult,Result,Model.Mesh);
Model.Set=Resume_set_information(Nset,Sset,Eset,Model.Mesh);
end