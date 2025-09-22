%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Construct_intersection_element(Model,id,varargin)
rangeX=[-inf,inf];rangeY=[-inf,inf];rangeZ=[-inf,inf];
%==========================Check input=====================================
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Sset'
            setid=varargin{vi+1};
        case 'RangeX'
            rangeX=varargin{vi+1};
        case 'RangeY'
            rangeY=varargin{vi+1};
        case 'RangeZ'
            rangeZ=varargin{vi+1};
        case 'Offset'
            offset=varargin{vi+1};
        case 'Meshdim'
            Meshdim=varargin{vi+1};
        otherwise
            disp(['Unknow input type is ignored:',varargin{vi}])
    end
end
%========================main function=====================================
nid=unique(Model.Mesh.surfaces(Model.Set.surface_set{setid(1)},:));
nid=nid(Model.Mesh.nodes(nid,1)>rangeX(1)&Model.Mesh.nodes(nid,1)<rangeX(2));
nid=nid(Model.Mesh.nodes(nid,2)>rangeY(1)&Model.Mesh.nodes(nid,2)<rangeY(2));
nid=nid(Model.Mesh.nodes(nid,3)>rangeZ(1)&Model.Mesh.nodes(nid,3)<rangeZ(2));
Model.Mesh.nodes(nid,:)=Model.Mesh.nodes(nid,:)+offset;
intersurface=[];
for si=2:1:numel(setid)
    nid2=unique(Model.Mesh.surfaces(Model.Set.surface_set{setid(si)},:));
    nid2=nid2(Model.Mesh.nodes(nid2,1)>rangeX(1)&Model.Mesh.nodes(nid2,1)<rangeX(2));
    nid2=nid2(Model.Mesh.nodes(nid2,2)>rangeY(1)&Model.Mesh.nodes(nid2,2)<rangeY(2));
    nid2=nid2(Model.Mesh.nodes(nid2,3)>rangeZ(1)&Model.Mesh.nodes(nid2,3)<rangeZ(2));
    cross_node=[nid;nid2];
    DT = delaunayTriangulation(Model.Mesh.nodes(cross_node,Meshdim));
    C=DT.ConnectivityList;
    intersurface=[intersurface;cross_node(C)];
end
Model.Mesh.surfaces=[Model.Mesh.surfaces;intersurface];
Model.Set.surface_set{id}=size(Model.Mesh.surfaces,1)-size(intersurface,1)+1:1:size(Model.Mesh.surfaces,1);
end