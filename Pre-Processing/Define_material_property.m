%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Define_material_property(Model,MID,varargin)
Key=Keywords();
%==========================Check input=====================================
set=zeros(0,1);Pset=zeros(0,1);
elasticity=struct('type',Key.linear,'parameter',zeros(0,1));
section=struct('type',Keywords().none,'parameter',zeros(0,1),'Update',false);%plate: H; %beam: A
plasticity=struct('type',Key.none,'parameter',zeros(0,1),'harden_parameter',zeros(0,3));
fluid=struct('type',Key.none,'conductivity',zeros(0,1),'weight',[0;3],'stabilization',10);
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Eset'
            set=varargin{vi+1};
        case 'Pset'
            Pset=varargin{vi+1};
        case 'Elasticity'
            elasticity=varargin{vi+1};
            elasticity.parameter=elasticity.parameter(:);
        case 'Section'
            section=varargin{vi+1};
            section.parameter=section.parameter(:);
        case 'Plasticity'
            plasticity=varargin{vi+1};
            plasticity.parameter=plasticity.parameter(:);
        case 'Fluid'
            fluid=varargin{vi+1};
            fluid.conductivity=fluid.conductivity(:);
            fluid.weight=fluid.weight(:);
        otherwise
            warning(['Unknow input type is ignored:',varargin{vi}])
    end
end
%========================main function=====================================
elasticity.D_matrix=zeros(0,0);
fluid.C_matrix=zeros(0,0);
%--------------------------------------------------------------------------
if size(plasticity.harden_parameter,1)>1
    H=plasticity.harden_parameter(2:end,:)-plasticity.harden_parameter(1:end-1,:);
    H=H(:,2)./H(:,1);H(end+1)=0;
    plasticity.harden_parameter(:,3)=H;
else
    plasticity.harden_parameter(:,3)=0;
end
%--------------------------------------------------------------------------
Model.Material_property(MID,1).set=set(:);
Model.Material_property(MID,1).particle_set=Pset(:);
Model.Material_property(MID,1).elasticity=elasticity;
Model.Material_property(MID,1).section=section;
Model.Material_property(MID,1).plasticity=plasticity;
Model.Material_property(MID,1).fluid=fluid;
end