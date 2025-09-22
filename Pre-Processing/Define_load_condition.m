%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function obj=Define_load_condition(obj,LID,varargin)
%==========================Check input=====================================
type=char(zeros(1,0));set=zeros(0,1);direction=zeros(0,1);
local_system=eye(size(obj.Mesh.nodes,2));expression=cell(0,1);
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
        case 'Direction'
            direction=varargin{vi+1};
        case 'Local_system'
            local_system=varargin{vi+1};
        otherwise
            warning(['Unknow input type is ignored:',varargin{vi}])
    end
end
if numel(expression)~=numel(direction)
    error('The size of Expression must be same as Direction');
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
obj.Load_condition(LID,1).type=type;
obj.Load_condition(LID,1).set=set(:);
obj.Load_condition(LID,1).value=value(:);
obj.Load_condition(LID,1).expression=expression(:);
obj.Load_condition(LID,1).direction=direction(:);
obj.Load_condition(LID,1).local_system=local_system;
end