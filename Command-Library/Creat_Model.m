%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Creat_Model()
Mesh=Creat_new_mesh();
%--------------------------------------------------------------------------
Set=struct('node_set',{cell(0,1)},'element_set',{cell(0,1)},'surface_set',{cell(0,1)});
%--------------------------------------------------------------------------
Material_property=struct('set',cell(0,1),'particle_set',cell(0,1),'elasticity',cell(0,1),'section',cell(0,1),'plasticity',cell(0,1));
%--------------------------------------------------------------------------
Boundary_condition=struct('type',cell(0,1),'set',cell(0,1),'value',cell(0,1),'expression',cell(0,1),'direction',cell(0,1));
Load_condition=struct('type',cell(0,1),'set',cell(0,1),'value',cell(0,1),'expression',cell(0,1),'direction',cell(0,1),'local_system',cell(0,1));
Contact_condition=struct('type',cell(0,1),'slave_set',cell(0,1),'master_set',cell(0,1),'constraint',cell(0,1),'discretization',cell(0,1),'enforcement',cell(0,1),'parameter',cell(0,1),'large_sliding',cell(0,1),'reseparation',cell(0,1));
Fluid_pressure_condition=struct('type',cell(0,1),'set',cell(0,1),'value',cell(0,1),'expression',cell(0,1));
Fluid_flux_condition=struct('type',cell(0,1),'set',cell(0,1),'value',cell(0,1),'expression',cell(0,1));
Predefined_field=struct('type',cell(0,1),'set',cell(0,1),'value',cell(0,1),'expression',cell(0,1),'used',cell(0,1));
%=============================Output=======================================
Model.Mesh=Mesh;
Model.Set=Set;
Model.Material_property=Material_property;
Model.Boundary_condition=Boundary_condition;
Model.Load_condition=Load_condition;
Model.Contact_condition=Contact_condition;
Model.Fluid_pressure_condition=Fluid_pressure_condition;
Model.Fluid_flux_condition=Fluid_flux_condition;
Model.Predefined_field=Predefined_field;
end