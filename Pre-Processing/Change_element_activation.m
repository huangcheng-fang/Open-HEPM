%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Change_element_activation(Model,varargin)
%==========================Check input=====================================
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Eset'
            eset=varargin{vi+1};
        case 'Status'
            status=varargin{vi+1};
        otherwise
            disp(['Unknow input type is ignored:',varargin{vi}])
    end
end
%========================main function=====================================
eid=Merge_cell(Model.Set.element_set,eset);
Model.Mesh.eactivation(eid)=logical(status);
nid=reshape(Model.Mesh.elements(Model.Mesh.eactivation,:),1,[]);
nid(isnan(nid))=[];
Model.Mesh.nactivation(:)=false;
Model.Mesh.nactivation(nid)=true;
%--------------------------------------------------------------------------
Model.Mesh=Get_surface(Model.Mesh);
end