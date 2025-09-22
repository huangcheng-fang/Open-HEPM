%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Define_contact_condition(Model,CID,varargin)
%==========================Check input=====================================
slave_set=zeros(0,1);master_set=zeros(0,1);plane_orientation=1;
type='solid&solid';discretization='surface_to_surface';enforcement='dual_lagrange';constraint='friction';
large_sliding=true;reseparation=true;Gap_measurement='configuration';fluid_exchange=false;perturbation=false;
parameter=struct('friction_coeff',0,'cohesion',0,'Emodulus',0,'Gmodulus',0,'permeability',0,'PDASS',0);
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Type'
            type=varargin{vi+1};
        case 'Slave_set'
            slave_set=varargin{vi+1};
        case 'Master_set'
            master_set=varargin{vi+1};
        case 'Plane_orientation'
            plane_orientation=varargin{vi+1};
        case 'Discretization'
            discretization=varargin{vi+1};
        case 'Enforcement'
            enforcement=varargin{vi+1};
        case 'Constraint'
            constraint=varargin{vi+1};
        case 'Parameter'
            parameter=varargin{vi+1};
        case 'Large_sliding'
            large_sliding=varargin{vi+1};
        case 'Reseparation'
            reseparation=varargin{vi+1};
        case 'Gap_measurement'
            Gap_measurement=varargin{vi+1};
        case 'Fluid_exchange'
            fluid_exchange=varargin{vi+1};
        case 'Perturbation'
            perturbation=varargin{vi+1};
        otherwise
            disp(['Unknow input type is ignored:',varargin{vi}])
    end
end
%--------------------------------------------------------------------------
if plane_orientation~=1&&plane_orientation~=-1
   error('plane orientation must be 1 or -1')
end
%========================main function=====================================
Model.Contact_condition(CID,1).type=type;
Model.Contact_condition(CID,1).slave_set=slave_set(:);
Model.Contact_condition(CID,1).master_set=master_set(:);
Model.Contact_condition(CID,1).plane_orientation=plane_orientation;
Model.Contact_condition(CID,1).discretization=discretization;
Model.Contact_condition(CID,1).enforcement=enforcement;
Model.Contact_condition(CID,1).constraint=constraint;
Model.Contact_condition(CID,1).parameter=parameter;
Model.Contact_condition(CID,1).large_sliding=large_sliding;
Model.Contact_condition(CID,1).reseparation=reseparation;
Model.Contact_condition(CID,1).Gap_measurement=Gap_measurement;
Model.Contact_condition(CID,1).fluid_exchange=fluid_exchange;
Model.Contact_condition(CID,1).perturbation=perturbation;
end