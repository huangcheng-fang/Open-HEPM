function [Result,Model]=Submit_To_FEM_Solver_Implicit(Result,Model,varargin)
if numel(varargin)~=2&&~isempty(varargin);error('Incorrect input');end
if numel(varargin)==2;T=varargin{2};else;T=0;end;Result.TimeStep=T;
%========================Solver control parameter==========================
GNL=Result.Solver_control.geometric_nonlinear;
Eqsolver=Result.Solver_control.equation_solver;
Tolerance=Result.Solver_control.Tolerance;
field_dof_index=[1,2];
%===========================Model parameter================================
Set=Model.Set;Material=Model.Material_property;
Boundary_condition=Model.Boundary_condition;
Load_condition=Model.Load_condition;
CD=Model.Contact_condition;
Predefined_field=Model.Predefined_field;
Mesh=Model.Mesh;
%===========================Field variable=================================
Nfield=Result.Nfield;
Efield=Result.Efield;
%============================Initiation====================================
[DOF,dof_range]=Get_node_dof_list(size(Mesh.nodes,1),{3,4});
max_dof=dof_range(2,1);
stepX=zeros(max_dof,1);
%--------------------------------------------------------------------------
MID=Get_element_material_ID(Material,Set.element_set,size(Mesh.elements,1));
uconstraint=Get_displacement_constraint(Boundary_condition,Set.node_set,Mesh.nodes,DOF{1});
% [Result,Predefined_field]=Form_predefined_field(Result,Predefined_field,Mesh,Set);
%--------------------------------------------------------------------------
[CP,CE]=Initial_contact_solver(numel(CD),Nfield.contact_stress,Nfield.contact_status);
CP=Update_contact_pair(CP,CD,Mesh,Set,stepX);
[CP,CE]=Execute_contact_analysis_new(CP,CE,CD,Mesh,Set,stepX,false);
%--------------------------------------------------------------------------
FEcell=Creat_FEcell(Mesh.etype,MID,1e5);
FEcell=Calculate_FEcell(FEcell,Mesh,field_dof_index);
FEdof=Assign_dof(FEcell.Elements,FEcell.Field_dof_num,DOF);
Efieldcell0=Assign_field_result(Efield,FEcell.Element_id);global Fint
%--------------------------------------------------------------------------
[Cdstrain,Cspin]=Get_element_dstrain(stepX,FEcell,FEdof(:,1),GNL);
[Efieldcell,Dmatrix]=Perform_constitutive_integration(Cdstrain,Cspin,Efieldcell0,FEcell.Material_id,Material,GNL);
Fint=Update_internal_force(Efieldcell.Stress,FEcell,FEdof(:,1),max_dof);
Fext=Form_F_global(Load_condition,Mesh,Set,DOF{1},max_dof);
%=========================Newton Iteration=================================
iternum=1;
while true
    Fun=Fext-Fint;
    K=Form_FEM_stiffness_matrix(FEcell,FEdof(:,1),Dmatrix,max_dof);
    [KC,FC]=Get_constrained_sysytem_new(K,Fun,uconstraint,Nfield.U,CE,DOF{1}(:,~Mesh.nactivation));
    %----------------------------------------------------------------------
    dX0=Equation_solver(KC,FC,Tolerance.PCGtol,Eqsolver);
    dX=CE.Dual_Lagrange.matrix*dX0;stepX=stepX+dX;
    Nfield.U=Nfield.U+dX;
    Mesh.positions=Mesh.positions+reshape(dX,3,[])';
    %----------------------------------------------------------------------
    [CP,CE]=Execute_contact_analysis_new(CP,CE,CD,Mesh,Set,stepX,false);
    %----------------------------------------------------------------------
    if GNL
        Mesh.nodes=Mesh.positions;
        Mesh=Update_mesh_geometry(Mesh);
        FEcell=Calculate_FEcell(FEcell,Mesh,field_dof_index);
        Fext=Form_F_global(Load_condition,Mesh,Set,DOF{1},max_dof);
    end
    [Cdstrain,Cspin]=Get_element_dstrain(stepX,FEcell,FEdof(:,1),GNL);
    [Efieldcell,Dmatrix]=Perform_constitutive_integration(Cdstrain,Cspin,Efieldcell0,FEcell.Material_id,Material,GNL);
    %----------------------------------------------------------------------
    Fint=Update_internal_force(Efieldcell.Stress,FEcell,FEdof(:,1),max_dof);
    max(abs(Nfield.U(2:3:end)))
    %----------------------------------------------------------------------
    Efield=Return_field_result(Efield,Efieldcell);
    Result.Nfield=Nfield;Result.Efield=Efield;
%     Plot_result(Result,Mesh,23,'S_mises',1,[]);drawnow
    Cvg=Is_Convergence_new(dX,stepX,Fext,Fint,uconstraint,CE,iternum,Tolerance);
    if Cvg;break;end;iternum=iternum+1;
end
%==============================Output======================================
Model.Mesh=Mesh;
Model.Predefined_field=Predefined_field;
Efield=Return_field_result(Efield,Efieldcell);
Result.Nfield=Nfield;Result.Efield=Efield;
% Result.NCflux=Result.NCflux/dT;
end