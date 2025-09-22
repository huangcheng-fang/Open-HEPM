%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Define_interface_element(Model,ID,varargin)
%==========================Check input=====================================
set=zeros(0,1);
parameter=struct('elastic_modulus',0,'poisson_ratio',0,'permeability',0,'flow_type',Keywords().interface_cubic);
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Sset'
            set=varargin{vi+1};
        case 'Parameter'
            parameter=varargin{vi+1};
        otherwise
            disp(['Unknow input type is ignored:',varargin{vi}])
    end
end
%========================main function=====================================
Model.Interface_element(ID,1).set=set(:);
Model.Interface_element(ID,1).parameter=parameter;
end