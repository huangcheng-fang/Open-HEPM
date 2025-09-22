%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Result,Model]=Submit_To_HEPM_Solver_Implicit(Result,Model,varargin)
if numel(varargin)~=2&&~isempty(varargin);error('Incorrect input');end
if numel(varargin)==2;T=varargin{2};else;T=0;end;Result.TimeStep=T;
fprintf("\n========================New Time Step========================\n\n");
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
[Material,Set.element_set,Set.particle_set,Part_list]=Split_material_by_part(Material,Set.element_set,Set.particle_set,Mesh.einpart,Mesh.pinpart);
%===========================Field variable=================================
Nfield=Result.Nfield;
Efield=Result.Efield;
Pfield0=Result.Pfield;
%============================Initiation====================================
[DOF,dof_range]=Get_node_dof_list(size(Mesh.nodes,1),{3});
max_dof=dof_range(end);
stepU=zeros(max_dof,1);
%--------------------------------------------------------------------------
EMID=Get_element_material_ID(Material,Set.element_set,size(Mesh.elements,1));
PMID=Get_particle_material_ID(Material,Set.particle_set,size(Mesh.particles,1));
uconstraint=Get_displacement_constraint(Boundary_condition,Set.node_set,Mesh.nodes,DOF{1});
%--------------------------------------------------------------------------
[CP,CE]=Initial_contact_solver(numel(CD),Nfield.contact_stress,Nfield.contact_status);
CP=Update_contact_pair(CP,CD,Mesh,Set,stepU);
[CP,CE]=Execute_contact_analysis_new(CP,CE,CD,Mesh,Set,stepU,false);
%--------------------------------------------------------------------------
FEcell=Creat_FEcell(Mesh.etype,EMID,1e50);
FEcell=Calculate_FEcell(FEcell,Mesh,field_dof_index);
FEdof=Assign_dof(FEcell.Elements,FEcell.Field_dof_num,DOF);
%--------------------------------------------------------------------------
[FPcell,CPfield0]=Creat_FPcell(PMID,EMID,Mesh.particles,Pfield0);
[FPcell,CPfield0]=Calculate_interaction_matrix(FPcell,CPfield0,FEcell);
FPcell=Update_particle_volume(FPcell,FEcell);
%--------------------------------------------------------------------------
[CEdstrain,CEspin]=Get_element_dstrain(stepU,FEcell,FEdof(:,1),GNL);
[CPdstrain,CPspin]=Get_particle_dstrain(CEdstrain,CEspin,FEcell.Material_id,FPcell,GNL);
[CPfield,Dmatrix]=Perform_constitutive_integration(CPdstrain,CPspin,CPfield0,FPcell.Material_id,Material,GNL);
CEStress=Get_gauss_point_stress(CPfield.Stress,FPcell,FEcell.DetJac,FEcell.Material_id);
%--------------------------------------------------------------------------
Fint=Update_internal_force(CEStress,FEcell,FEdof(:,1),max_dof);
Fext=Form_F_global(Load_condition,Mesh,Set,DOF{1},max_dof);
% C=Form_stiffness_matrix_FPEM(FEcell,FEdof(:,1),FPcell,Dmatrix,max_dof)*0.8;
%=========================Newton Iteration=================================
iternum=1;
while true
    fprintf("\n------------------Newton Iteration %d------------------\n\n",iternum);
    Fun=Fext-Fint;
    K=Form_stiffness_matrix_FPEM(FEcell,FEdof(:,1),FPcell,Dmatrix,max_dof);
%     K=K+C;
    [KC,FC]=Get_constrained_sysytem_new(K,Fun,uconstraint,Nfield.U,CE,DOF{1}(:,~Mesh.nactivation));
    %----------------------------------------------------------------------
    dU0=Equation_solver(KC,FC,Tolerance.PCGtol,Eqsolver);
    dU=CE.Dual_Lagrange.matrix*dU0;
    stepU=stepU+dU;Nfield.U=Nfield.U+dU;
    Mesh.positions=Mesh.positions+reshape(dU,3,[])';
    %----------------------------------------------------------------------
    [CP,CE]=Execute_contact_analysis_new(CP,CE,CD,Mesh,Set,stepU,iternum>20);
    %----------------------------------------------------------------------
    if GNL
        Mesh.nodes=Mesh.positions;
        Mesh=Update_mesh_geometry(Mesh);
        FPcell=Update_particle_position(FPcell,FEcell,DOF{1},dU);
        FEcell=Calculate_FEcell(FEcell,Mesh,field_dof_index);
    end
    %--------------------------------------------------------------------------  
    FPcell=Update_particle_volume(FPcell,FEcell);
    %--------------------------------------------------------------------------
    [CEdstrain,CEspin]=Get_element_dstrain(stepU,FEcell,FEdof(:,1),GNL);
    [CPdstrain,CPspin]=Get_particle_dstrain(CEdstrain,CEspin,FEcell.Material_id,FPcell,GNL);
    [CPfield,Dmatrix]=Perform_constitutive_integration(CPdstrain,CPspin,CPfield0,FPcell.Material_id,Material,GNL);
    CEStress=Get_gauss_point_stress(CPfield.Stress,FPcell,FEcell.DetJac,FEcell.Material_id);
    %----------------------------------------------------------------------
    Fint=Update_internal_force(CEStress,FEcell,FEdof(:,1),max_dof);
%     Fint=Fint+C*stepU;
    %----------------------------------------------------------------------
    Efieldcell=Convert_variable_FPEM2FEM(CPfield,FPcell,FEcell.DetJac,FEcell.Material_id,FEcell.Element_id);
    Efield=Return_field_result(Efield,Efieldcell);
    Result.Nfield=Nfield;Result.Efield=Efield;
    Plot_result(Result,Mesh,23,'Sy',1,[]);drawnow
    %----------------------------------------------------------------------
    Cvg=Is_Convergence_new(dU,stepU,Fext,Fint,uconstraint,CE,iternum,Tolerance);
    if Cvg;break;end;iternum=iternum+1;
end
Efieldcell=Convert_variable_FPEM2FEM(CPfield,FPcell,FEcell.DetJac,FEcell.Material_id,FEcell.Element_id);
Efield=Return_field_result(Efield,Efieldcell);
[Mesh.particles,Mesh.pinpart]=Return_particle_information(FPcell,Part_list);
Pfield=Return_field_result_FPEM(CPfield);
%==============================Output======================================
Model.Mesh=Mesh;
Model.Predefined_field=Predefined_field;
Result.Nfield=Nfield;
Result.Efield=Efield;
Result.Pfield=Pfield;
end