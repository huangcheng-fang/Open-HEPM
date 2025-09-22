%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Geometry=Generate_geometry(Mesh,partId,addNodes)
nodes=Mesh.nodes;
surfaces=Mesh.surfaces;
Sinpart=Mesh.sinpart;
Stype=Mesh.stype;
if size(surfaces,2)<10;surfaces(:,10)=0;end
%==========================================================================
loc=Sinpart==partId;
part_surfaces=surfaces(loc,:);
part_stype=Stype(loc);
loc3=part_stype==3;
loc4=part_stype==4;
part_surfaces=[part_surfaces(loc3,1:3);
               part_surfaces(loc4,1:3);
               part_surfaces(loc4,[3,4,1])];
Geometry=geometryFromMesh(createpde(),nodes',part_surfaces');
if ~isempty(addNodes)
    addVertex(Geometry,"Coordinates",addNodes);
end
% pdegplot(Geometry,"EdgeLabels","on","FaceLabels","on")
end