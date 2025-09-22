%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Result,Model]=Submit_To_USFEM_Solver(Result,Model,varargin)
if numel(varargin)~=2&&~isempty(varargin);error('Incorrect input');end
if numel(varargin)==2;T=varargin{2};else;T=0;end
dT=T-Result.Time_step;Result.Time_step=T;
%==============================Input=======================================
Set=Model.Set;Material=Model.Material_property;
Boundary_condition=Model.Boundary_condition;
Load_condition=Model.Load_condition;
CD=Model.Contact_condition;
Predefined_field=Model.Predefined_field;
Mesh=Model.Mesh;
maxdof=numel(Mesh.nodes);
GNL=Result.GeometricNonlinear;
EQSolver=Result.EquationSolver;
Tolerance=Result.Tolerance;
%============================Initiation====================================
[Material,MID]=Initialize_material(Material,Set.element_set,Mesh.eactivation);
[Result,Predefined_field]=Form_predefined_field(Result,Predefined_field,Mesh,Set);
Result.dU=zeros(size(Result.U));
Result.stepU=zeros(size(Result.U,1),1);
NCstress=Result.NCstress;
NCstatus=Result.NCstatus;
Result.Estress0=Result.Estress;
Result.Estrain0=Result.Estrain;
Result.Epstrain0=Result.Epstrain;
Result.Eparameter0=Result.Eparameter;
%--------------------------------------------------------------------------
FE=FEM_Element(Mesh,MID);
FE=Assign_section_parameter(FE,Material);
[Result,FE]=Update_element_field(Result,Mesh,FE,Material,GNL);
%--------------------------------------------------------------------------
[CP,CE]=Initial_contact_solver(numel(CD),NCstress,NCstatus);
CP=Update_contact_pair(CP,CD,Mesh,Set,Result.stepU);
[CP,CE]=Execute_contact_analysis_new(CP,CE,CD,Mesh,Set,Result.stepU,false);
%--------------------------------------------------------------------------
Result.FEinfo=struct('type',{FE.type}','eid',{FE.eid}','detJac',{FE.detJac}','dof',{FE.dof}','N_matrix',{FE.N_matrix}');
NResult=Mapping_variable_E2N(Result);
Resulttemp=Mapping_variable_N2E(NResult,Result,Mesh);
Fint=Update_internal_force(Resulttemp,FE,maxdof);
Fext=Form_F_global(Load_condition,Mesh,Set,maxdof);
%--------------------------------------------------------------------------
Uconstraint=Get_displacement_constraint(Boundary_condition,Set.node_set,Mesh.nodes,Result.U);
%=========================Newton Iteration=================================
iternum=1;
while 1
    Fun=Fext-Fint;
    K_semi=Form_USFEM_stiffness_matrix(FE,maxdof);whos("K_semi")
    %----------------------------------------------------------------------
    tic;[X,iter]=JPCG_for_USFEM(K_semi,Fun,[],Tolerance.PCGtol,10000,[],[],Uconstraint);toc,iter
    Result=Split_result_new(Result,CE.Dual_Lagrange.matrix*X,maxdof);
    Mesh.positions=Mesh.positions+reshape(Result.dU,3,[]).';
    if iternum==1;Uconstraint(:,2)=0;end
    %----------------------------------------------------------------------
    CP=Get_contact_stress_new(CP,CD,Fun-K_semi*(K_semi'*Result.dU),[],[]);
    [CP,CE]=Execute_contact_analysis_new(CP,CE,CD,Mesh,Set,Result.stepU,false);
    %----------------------------------------------------------------------
    if GNL
        Mesh.nodes=Mesh.nodes+reshape(Result.dU,3,[]).';
        Mesh=Update_mesh_geometry(Mesh);
        FE=FEM_Element(Mesh,MID);
        Fext=Form_F_global(Load_condition,Mesh,Set,maxdof);
    end
    [Result,FE]=Update_element_field(Result,Mesh,FE,Material,GNL);
    %----------------------------------------------------------------------
    Result.FEinfo=struct('type',{FE.type}','eid',{FE.eid}','detJac',{FE.detJac}','dof',{FE.dof}','N_matrix',{FE.N_matrix}');
    NResult=Mapping_variable_E2N(Result);
    Resulttemp=Mapping_variable_N2E(NResult,Result,Mesh);
    Fint=Update_internal_force(Resulttemp,FE,maxdof);
    %----------------------------------------------------------------------
    try close(23);catch;end
    Plot_result(Result,Mesh,23,'Uz',2,[]);drawnow
    %----------------------------------------------------------------------
    Cvg=Is_Convergence_new(Result,Fext,Fint,Uconstraint,CE,iternum,Tolerance);
    if Cvg;break;end;iternum=iternum+1;
end
Model.Mesh=Mesh;
Model.Predefined_field=Predefined_field;
Result.FEinfo=struct('type',{FE.type}','eid',{FE.eid}','detJac',{FE.detJac}','dof',{FE.dof}','N_matrix',{FE.N_matrix}');
end