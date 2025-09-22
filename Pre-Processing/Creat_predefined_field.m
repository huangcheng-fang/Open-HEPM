%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Creat_predefined_field(Model,ID,varargin)
%==========================Check input=====================================
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Type'
            type=varargin{vi+1};
        case 'Eset'
            set=varargin{vi+1};
        case 'Nset'
            set=varargin{vi+1};
        case 'Expression'
            expression=varargin{vi+1};
        otherwise
            warning(['Unknow input type is ignored:',varargin{vi}])
    end
end
%--------------------------------------------------------------------------
switch type
    case 'stress'
        if numel(expression)~=6
            error('Predefined stress must input six expressions:{Sx,Sy,Sz,Sxy,Syz,Szx}')
        end
    case 'void_ratio'
        if numel(expression)~=1
            error('Predefined void ratio can only input one expression')
        end
    case 'consolidation_pressure'
        if numel(expression)~=1
            error('Predefined consolidation pressure can only input one expression')
        end
    otherwise
        error('Unknown field type')
end
%--------------------------------------------------------------------------
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
Model.Predefined_field(ID,1).type=type;
Model.Predefined_field(ID,1).set=set(:);
Model.Predefined_field(ID,1).value=value(:);
Model.Predefined_field(ID,1).expression=expression(:);
Model.Predefined_field(ID,1).used=false;
end