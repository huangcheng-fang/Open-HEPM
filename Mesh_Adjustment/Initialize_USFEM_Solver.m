%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Result=Initialize_USFEM_Solver(Model,varargin)
%==========================Check input=====================================
Parallelization=1;GeometricNonlinear=true;EquationSolver='default';
Tolerance=struct('Utol',1e-5,'Ftol',1e-5,'PCGtol',1e-8,'MaxIterNum',100);
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Parallelization'
            Parallelization=varargin{vi+1};
        case 'GeometricNonlinear'
            GeometricNonlinear=varargin{vi+1};
        case 'Tolerance'
            Tolerance=varargin{vi+1};
        case 'EquationSolver'
            EquationSolver=varargin{vi+1};
        otherwise
            warning(['Unknow input type is ignored:',varargin{vi}])
    end
end
%========================main function=====================================
node_num=size(Model.Mesh.nodes,1);
element_num=size(Model.Mesh.elements,1);
dimension=size(Model.Mesh.nodes,2);
maxdof=node_num*dimension;
EIPnum=8;

Result.U=zeros(maxdof,1);
Result.dU=zeros(maxdof,1);
Result.stepU=zeros(maxdof,1);
Result.P=zeros(node_num,1);
Result.dP=zeros(node_num,1);

Result.Nstress=zeros(dimension*2,node_num);
Result.Nstrain=zeros(dimension*2,node_num);
Result.Npstrain=zeros(dimension*2,node_num);
Result.Ndvstrain=zeros(dimension*2,node_num);
Result.Nparameter=nan(dimension*0,dimension*0,node_num);


Result.NCstress=zeros(maxdof,numel(Model.Contact_condition));
Result.NCflux=zeros(node_num,numel(Model.Contact_condition));

Result.Fext=zeros(maxdof,1);
Result.Time_step=0;
%--------------------------------------------------------------------------
Result.Tolerance=Tolerance;
Result.EquationSolver=EquationSolver;
Result.GeometricNonlinear=GeometricNonlinear;
Result.ParallelWorkerNum=Start_Parallel_Processor(Parallelization);
Result.FEinfo=struct('type',cell(0,1),'eid',cell(0,1),'detJac',cell(0,1),'dof',cell(0,1),'N_matrix',cell(0,1));
end