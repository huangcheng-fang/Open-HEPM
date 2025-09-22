%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Result,Model]=Submit_To_SPFEM_Solver_Implicit(Result,Model,varargin)
if numel(varargin)~=2&&~isempty(varargin);error('Incorrect input');end
if numel(varargin)==2;T=varargin{2};else;T=0;end;Result.TimeStep=T;
%========================Solver control parameter==========================
GNL=Result.Solver_control.geometric_nonlinear;
Eqsolver=Result.Solver_control.equation_solver;
Tolerance=Result.Solver_control.Tolerance;
field_dof_index=[1];
%===========================Model parameter================================
Set=Model.Set;Material=Model.Material_property;
Boundary_condition=Model.Boundary_condition;
Load_condition=Model.Load_condition;
CD=Model.Contact_condition;
Predefined_field=Model.Predefined_field;
Mesh=Model.Mesh;
global Fint
%===========================Field variable=================================
Nfield=Result.Nfield;
Pfield0=Result.Pfield;
%============================Initiation====================================
[DOF,dof_range]=Get_node_dof_list(size(Mesh.nodes,1),{3});
max_dof=dof_range(end);
stepU=zeros(max_dof,1);
%--------------------------------------------------------------------------
MID=Get_particle_material_ID_SPFEM(Material,Set.element_set,Mesh);
uconstraint=Get_displacement_constraint(Boundary_condition,Set.node_set,Mesh.nodes,DOF{1});
% [Result,Predefined_field]=Form_predefined_field(Result,Predefined_field,Mesh,Set);
%--------------------------------------------------------------------------
[CP,CE]=Initial_contact_solver(numel(CD),Nfield.contact_stress,Nfield.contact_status);
CP=Update_contact_pair(CP,CD,Mesh,Set,stepU);
[CP,CE]=Execute_contact_analysis_new(CP,CE,CD,Mesh,Set,stepU,false);
%--------------------------------------------------------------------------
FEcell=FEM_Element(Mesh,field_dof_index);
FEdof=Assign_dof(FEcell.Elements,FEcell.Field_dof_num,DOF);
%--------------------------------------------------------------------------
[CEdstrain,CEspin]=Get_element_dstrain(stepU,FEcell,FEdof(:,1),GNL);
[Fai,node_volume]=Calculate_smoothing_factor(Mesh,CEdstrain,FEcell.Element_id);
[Pfield,Dmatrix]=Update_particle_field(CEdstrain,CEspin,Pfield0,Fai,FEcell,Material,MID,GNL);
Fint=Update_internal_force_SPFEM(Pfield.Stress,Fai,node_volume,FEcell,FEdof(:,1),max_dof);
Fext=Form_F_global(Load_condition,Mesh,Set,DOF{1},max_dof);
%=========================Newton Iteration=================================
iternum=1;
while true
    K=Form_SPFEM_stiffness_matrix(FEcell,FEdof(:,1),Dmatrix,Fai,node_volume,max_dof);
    Fun=Fext-Fint;
    [KC,FC]=Get_constrained_sysytem_new(K,Fun,uconstraint,Nfield.U,CE,DOF{1}(:,~Mesh.nactivation));
    %----------------------------------------------------------------------
    dU=KC\FC;
%     dU=Equation_solver(KC,FC,Tolerance.PCGtol,Eqsolver);
    dU=CE.Dual_Lagrange.matrix*dU;
    stepU=stepU+dU;Nfield.U=Nfield.U+dU;
    Mesh.positions=Mesh.positions+reshape(dU,3,[])';
    %----------------------------------------------------------------------
    CP=Get_contact_stress_new(CP,CD,Fun-K*dU,[],[]);
    [CP,CE]=Execute_contact_analysis_new(CP,CE,CD,Mesh,Set,stepU,false);
    %----------------------------------------------------------------------
    if GNL
        Mesh.nodes=Mesh.positions;
        Mesh=Update_mesh_geometry(Mesh);
        FEcell=FEM_Element(Mesh,field_dof_index);
        Fext=Form_F_global(Load_condition,Mesh,Set,DOF{1},max_dof);
    end
    [CEdstrain,CEspin]=Get_element_dstrain(stepU,FEcell,FEdof(:,1),GNL);
    [Fai,node_volume]=Calculate_smoothing_factor(Mesh,CEdstrain,FEcell.Element_id);
    [Pfield,Dmatrix]=Update_particle_field(CEdstrain,CEspin,Pfield0,Fai,FEcell,Material,MID,GNL);
    %----------------------------------------------------------------------
    Fint=Update_internal_force_SPFEM(Pfield.Stress,Fai,node_volume,FEcell,FEdof(:,1),max_dof);
    max(abs(Nfield.U(2:3:end)))
    %----------------------------------------------------------------------
    Result.Nfield=Nfield;Result.Pfield=Pfield;
%     Plot_result(Result,Mesh,23,'S_mises',1,[]);drawnow
    Cvg=Is_Convergence_new(dU,stepU,Fext,Fint,uconstraint,CE,iternum,Tolerance);
    if Cvg;break;end;iternum=iternum+1;
end
%==============================Output======================================
Model.Mesh=Mesh;
Model.Predefined_field=Predefined_field;
Result.Nfield=Nfield;
Result.Pfield=Pfield;
% Result.NCflux=Result.NCflux/dT;
end