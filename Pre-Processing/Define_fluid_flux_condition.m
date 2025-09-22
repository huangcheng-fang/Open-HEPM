%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Define_fluid_flux_condition(Model,FID,varargin)
%==========================Check input=====================================
set=zeros(0,1);expression=cell(0,1);
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Type'
            type=varargin{vi+1};
        case 'Nset'
            set=varargin{vi+1};
        case 'Eset'
            set=varargin{vi+1};
        case 'Sset'
            set=varargin{vi+1};
        case 'Expression'
            expression=varargin{vi+1};
        otherwise
            disp(['Unknow input type is ignored:',varargin{vi}])
    end
end
value=nan(numel(expression),1);
for ii=1:1:numel(expression)
    if isnumeric(expression{ii})
        value(ii)=expression{ii};
        expression{ii}=char(zeros(1,0));
    end
    if ~isnan(str2double(expression{ii}))
        value(ii)=str2double(expression{ii});
        expression{ii}=char(zeros(1,0));
    end
end
%========================main function=====================================
Model.Fluid_flux_condition(FID,1).type=type(:);
Model.Fluid_flux_condition(FID,1).set=set(:);
Model.Fluid_flux_condition(FID,1).value=value(:);
Model.Fluid_flux_condition(FID,1).expression=expression(:);
end