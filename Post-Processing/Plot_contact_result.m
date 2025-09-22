%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Plot_contact_result(Result,Model,n,type,scale,LineStyle,Permutation)
if nargin<6
    LineStyle='-';
end
if nargin<7
    Permutation=[1,2,3];
end
%--------------------------------------------------------------------------
surfaces=Model.Mesh.surfaces;
Stype=Model.Mesh.stype;
Sset=Model.Set.surface_set;
Contact_condition=Model.Contact_condition;
NCstress=Result.NCstress;
NCstressN=NCstress(3:3:end,:);
%--------------------------------------------------------------------------
switch type
    case 'CSn'
        value=NCstress(3:3:end,:);
    case 'CSt'
        value=sqrt(NCstress(1:3:end,:).^2+NCstress(2:3:end,:).^2);
    otherwise
        error(['Wrong Keyword: ',type]);
end
%==============================main function===============================
nodes=Model.Mesh.nodes+(scale-1)*reshape(Result.U,[],size(Model.Mesh.nodes,1)).';
nodes=nodes(:,Permutation);
figure(n)
for ci=1:1:numel(Contact_condition)
    slave=Merge_cell(Sset,Contact_condition(ci).slave_set);
    eid=slave(Stype(slave)==3);
    Tri=surfaces(eid,1:3);
    nodeID=unique(Tri);
    nodeID(isnan(nodeID))=[];
    nodeID=nodeID(NCstressN(nodeID,ci)~=0);
    Tri=Tri(all(ismember(Tri,nodeID),2),:);
    patch('Faces',Tri,'Vertices',nodes,'FaceVertexCData',value(:,ci),'FaceColor','interp','LineWidth',0.001,'LineStyle',LineStyle);
end
%--------------------------------------------------------------------------
axis equal
if size(nodes,2)==3
    view(-18,15)
end
colorbar
colormap(jet(18))
drawnow
end