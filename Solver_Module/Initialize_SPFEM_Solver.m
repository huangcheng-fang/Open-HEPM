%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Result=Initialize_SPFEM_Solver(Model,varargin)
%==========================Check input=====================================
parallelization=1;geometric_nonlinear=true;equation_solver='default';
Tolerance=struct('Utol',1e-5,'Ftol',1e-5,'PCGtol',1e-8,'MaxIterNum',100);
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Parallelization'
            parallelization=varargin{vi+1};
        case 'GeometricNonlinear'
            geometric_nonlinear=varargin{vi+1};
        case 'Tolerance'
            Tolerance=varargin{vi+1};
        case 'EquationSolver'
            equation_solver=varargin{vi+1};
        otherwise
            warning(['Unknow input type is ignored:',varargin{vi}])
    end
end
%========================main function=====================================
node_num=size(Model.Mesh.nodes,1);
dimension=size(Model.Mesh.nodes,2);
maxdof=node_num*dimension;
%--------------------------------------------------------------------------
Result.Time_step=0;
Result.Nfield.U=zeros(maxdof,1);
Result.Nfield.V=zeros(maxdof,1);
Result.Nfield.A=zeros(maxdof,1);
Result.Nfield.contact_stress=zeros(maxdof,numel(Model.Contact_condition));
Result.Nfield.contact_status=false(maxdof,numel(Model.Contact_condition));
Result.Nfield.contact_flux=zeros(node_num,numel(Model.Contact_condition));
%--------------------------------------------------------------------------
Result.Pfield.Stress=zeros(6,1,6,node_num);
Result.Pfield.Strain=zeros(6,1,6,node_num);
Result.Pfield.Pstrain=zeros(6,1,6,node_num);
Result.Pfield.Parameter=zeros(2,1,6,node_num);
%--------------------------------------------------------------------------
Result.Solver_control.Tolerance=Tolerance;
Result.Solver_control.equation_solver=equation_solver;
Result.Solver_control.geometric_nonlinear=geometric_nonlinear;
Result.Solver_control.parallel_worker_num=Start_Parallel_Processor(parallelization);
end