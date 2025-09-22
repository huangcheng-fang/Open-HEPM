%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Plot_surface(Mesh,n)
if isfield(Mesh,'Mesh')
    Mesh=Mesh.Mesh;
end
nodes=Mesh.nodes;
elements=Mesh.elements;
surfaces=Mesh.surfaces;
sinpart=Mesh.sinpart;
sinfacet=Mesh.sinfacet;
einpart=Mesh.einpart;
etype=Mesh.etype;
%==============================main function===============================
figure(n)
facet=unique(sinfacet);
faceColor=colorcube(numel(facet));
colorbar
colormap(faceColor)
for ii=1:1:size(facet,1)
    patch('Faces',surfaces(sinfacet==facet(ii),:),'Vertices',nodes,'FaceVertexCData',facet(ii),'FaceColor',faceColor(ii,:),'LineWidth',1.0);
end
% facet=unique(einpart);
% for ii=1:1:size(facet,1)
%     patch('Faces',elements((etype==223|etype==312)&einpart==facet(ii),1:end-1),'Vertices',nodes,'FaceColor',[38,255,255]/255,'LineWidth',1.0);
% end
%--------------------------------------------------------------------------
if size(nodes,2)==3
    view(-18,15)
end
axis equal
drawnow
end