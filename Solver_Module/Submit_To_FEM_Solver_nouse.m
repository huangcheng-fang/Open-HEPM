%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Result,Model]=Submit_To_FEM_Solver(Result,Model,varargin)
if numel(varargin)~=2&&~isempty(varargin);error('Incorrect input');end
if numel(varargin)==2;T=varargin{2};else;T=0;end
dT=T-Result.Time_step;Result.Time_step=T;
%==============================Input=======================================
Set=Model.Set;Material=Model.Material_property;
Boundary_condition=Model.Boundary_condition;
Load_condition=Model.Load_condition;
Contact_condition=Model.Contact_condition;
Fluid_pressure_condition=Model.Fluid_pressure_condition;
Fluid_flux_condition=Model.Fluid_flux_condition;
Predefined_field=Model.Predefined_field;
Mesh=Model.Mesh;
maxdof=numel(Mesh.nodes);
GNL=Result.GeometricNonlinear;
EQSolver=Result.EquationSolver;
Tolerance=Result.Tolerance;
%============================Initiation====================================
Result.dU=zeros(size(Result.dU));
Result.dP=zeros(size(Result.dP));
Result.stepU=zeros(size(Result.U,1),1);
Result.NCflux=Result.NCflux*dT;
iternum=1;Pprevious=Result.P;
[Material,Contact_condition]=Permeability_X_Time(Material,Contact_condition,dT);
[Material,MID]=Initialize_material(Material,Set.element_set,Mesh.eactivation);
[Result,Predefined_field]=Form_predefined_field(Result,Predefined_field,Mesh,Set);
%--------------------------------------------------------------------------
FE=FEM_Element(Mesh,MID);
FE=Assign_section_parameter(FE,Material);
[Result,FE]=Update_element_field(Result,Mesh,FE,Material,GNL);
[Result,FE]=Update_element_fluid_field(Result,FE,Material,GNL);
%--------------------------------------------------------------------------
Fint=Update_internal_force(Result,FE,maxdof);
Flux_couple=Update_Flux_couple(Result,FE,maxdof);
Fext=Form_F_global(Load_condition,Mesh,Set,maxdof);
Flux_ext=Form_fluid_flux_global(Fluid_flux_condition,Mesh,Set)*dT;
Fluid_gravity=Form_FEM_fluid_gravity_vector(FE,maxdof);
%--------------------------------------------------------------------------
CP=Initial_contact_solver(numel(Contact_condition));
CP=Update_contact_pair(CP,Contact_condition,Mesh,Set);
Result=Get_contact_stress(Result,CP,Contact_condition,Fext-Fint,[]);
[CP,CE]=Execute_contact_analysis(CP,Contact_condition,Mesh,Result,false,true);
%--------------------------------------------------------------------------
Uconstraint=Get_displacement_constraint(Boundary_condition,Set.node_set,Mesh.nodes,Result.U);
Pconstraint=Get_fluid_pressure_constraint(Fluid_pressure_condition,Set.node_set,Mesh.nodes,Result.P);
constraint=[Uconstraint;Pconstraint];
%=========================Newton Iteration=================================
while 1
    KS=Form_FEM_stiffness_matrix(FE,maxdof);
    [KF,KCP]=Form_FEM_fluid_matrix(FE,maxdof);
    PPP=PPP_stabilization(FE,Material,maxdof);
    Fint=Fint-KCP*Result.P;Flux_int=Flux_couple+KF*Result.P+PPP*(Result.P-Pprevious)+Fluid_gravity;
    KHM=[KS,-KCP;-KCP.',-KF-PPP];FHM=[Fext-Fint;Flux_int-Flux_ext];
    %----------------------------------------------------------------------
    [KHMC,FHMC]=Get_constrained_sysytem(KHM,FHM,constraint,CE);
    if iternum==1;constraint(:,2)=0;end
    %----------------------------------------------------------------------
    X=Equation_solver(KHMC,FHMC,Tolerance.PCGtol,EQSolver);
    Result=Split_result(Result,CE.Dual_Lagrange.matrix*X(1:maxdof+maxdof/3),maxdof);
    %----------------------------------------------------------------------
    if GNL
        Mesh.nodes=Mesh.nodes+reshape(Result.dU,3,[]).';
        Mesh=Update_mesh_geometry(Mesh);
        FE=FEM_Element(Mesh);
        Fext=Form_F_global(Load_condition,Mesh,Set,maxdof);
        Flux_ext=Form_fluid_flux_global(Fluid_flux_condition,Mesh,Set)*dT;
        Fluid_gravity=Form_FEM_fluid_gravity_vector(FE,maxdof);
    else
        Mesh.positions=Mesh.positions+reshape(Result.dU,3,[]).';
    end
    %----------------------------------------------------------------------
    if norm(Result.dU)/norm(Result.stepU)>1e-2;CP=Update_contact_pair(CP,Contact_condition,Mesh,Set);end
    Result=Get_contact_stress(Result,CP,Contact_condition,FHM-KHM*[Result.dU;Result.dP],X(maxdof+maxdof/3+1:end));
    [CP,CE]=Execute_contact_analysis(CP,Contact_condition,Mesh,Result,false,true);
    %----------------------------------------------------------------------
    [Result,FE]=Update_element_field(Result,Mesh,FE,Material,GNL);
    [Result,FE]=Update_element_fluid_field(Result,FE,Material,GNL);
    %----------------------------------------------------------------------
    Fint=Update_internal_force(Result,FE,maxdof);
    Flux_couple=Update_Flux_couple(Result,FE,maxdof);
    %----------------------------------------------------------------------
    try close(23);catch;end
    Result.FEinfo=struct('type',{FE.type}','eid',{FE.eid}','detJac',{FE.detJac}','dof',{FE.dof}','N_matrix',{FE.N_matrix}');
    Plot_result(Result,Mesh,23,'Uz',2,[]);drawnow
    %----------------------------------------------------------------------
    Fint_temp=Fint-KCP*Result.P;
    Flux_int_temp=Flux_couple+KF*Result.P+PPP*(Result.P-Pprevious)+Fluid_gravity;
    Cvg=Is_Convergence(Result,Fext,Fint_temp,Flux_ext,Flux_int_temp,constraint,CE,iternum,Tolerance);
    if Cvg;break;end;iternum=iternum+1;
end
Model.Mesh=Mesh;
Model.Predefined_field=Predefined_field;
Result.FEinfo=struct('type',{FE.type}','eid',{FE.eid}','detJac',{FE.detJac}','dof',{FE.dof}','N_matrix',{FE.N_matrix}');
Result.NCflux=Result.NCflux/dT;
end