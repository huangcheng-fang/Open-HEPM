%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Construct_virtual_surface(Model,varargin)
%==========================Check input=====================================
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Sset'
            setid=varargin{vi+1};
        otherwise
            disp(['Unknow input type is ignored:',varargin{vi}])
    end
end
%========================main function=====================================
sid=Merge_cell(Model.Set.surface_set,setid);
surface=Model.Mesh.surfaces(sid,:);
[nid,~,newsurface]=unique(surface);
newsurface=reshape(newsurface,[],size(surface,2))+size(Model.Mesh.nodes,1);
Model.Mesh.nodes=[Model.Mesh.nodes;Model.Mesh.nodes(nid,:)];
Model.Mesh.ninpart=[Model.Mesh.ninpart;repmat(max(Model.Mesh.ninpart)+1,numel(nid),1)];
Model.Mesh.surfaces=[Model.Mesh.surfaces;newsurface];
Model.Mesh.stype=[Model.Mesh.stype;Model.Mesh.stype(sid)];
Model.Mesh.sinpart=[Model.Mesh.sinpart;repmat(max(Model.Mesh.sinpart)+1,numel(sid),1)];
Model.Mesh.sinfacet=[Model.Mesh.sinfacet;repmat(max(Model.Mesh.sinfacet)+1,numel(sid),1)];
Model.Mesh.scenter=[Model.Mesh.scenter;Model.Mesh.scenter(sid,:)];
end