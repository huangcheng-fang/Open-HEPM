%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Import_Mesh(Model,varargin)
%==========================Check input=====================================
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'File'
            filename=varargin{vi+1};
        otherwise
            error(['Unknow input type:',varargin{vi}])
    end
end
%========================main function=====================================
Model.Mesh=Import_mesh(Model.Mesh,filename);
if isempty(Model.Mesh.surfaces);Model.Mesh=Get_surface(Model.Mesh);end
Model.Mesh=Update_mesh_geometry(Model.Mesh);
end