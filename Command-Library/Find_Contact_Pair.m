%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [pair,orient]=Find_Contact_Pair(Model,varargin)
%==========================Check input=====================================
slave_part=unique(Model.Mesh.sinpart);master_part=slave_part;plotfig='off';tolgap=0.1;
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'RangeSlavePart'
            slave_part=varargin{vi+1};
        case 'RangeMasterPart'
            master_part=varargin{vi+1};
        case 'RangePart'
            slave_part=varargin{vi+1};
            master_part=varargin{vi+1};
        case 'IsPlot'
            plotfig=varargin{vi+1};
        case 'Tolgap'
            tolgap=varargin{vi+1};
        otherwise
            error(['Wrong input argument:',command{1}]);
    end
end
%========================main function=====================================
Mesh=Model.Mesh;
nodes=Mesh.nodes;surfaces=Mesh.surfaces;Scenter=Mesh.scenter;Sinpart=Mesh.sinpart;Sinfacet=Mesh.sinfacet;
[pair,orient]=Find_contact_pair(slave_part(:),master_part(:),nodes,surfaces,Scenter,Sinpart,Sinfacet,tolgap);
figure;hold on
if plotfig
    for ii=1:1:size(pair,1)
        color=rand(1,3);nid=unique(surfaces(Sinfacet==pair(ii,1),:));
        patch('Faces',surfaces(Sinfacet==pair(ii,1),:),'Vertices',nodes,'FaceColor',color,'Facealpha',0.8);
        patch('Faces',surfaces(Sinfacet==pair(ii,2),:),'Vertices',nodes,'FaceColor',color,'Facealpha',0.8);
        plot3(nodes(nid,1),nodes(nid,2),nodes(nid,3),'.','MarkerSize',10)
    end
view(-18,15);drawnow;axis equal
end
end