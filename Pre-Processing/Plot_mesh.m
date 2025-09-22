%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Plot_mesh(Mesh,n)
nodes=Mesh.nodes;
elements=Mesh.elements;
surfaces=Mesh.surfaces;
sinpart=Mesh.sinpart;
stype=Mesh.stype;
surfaces(stype==6,4:end)=nan;
einpart=Mesh.einpart;
etype=Mesh.etype;
%==============================main function===============================
figure(n)
part=unique(sinpart);
for ii=1:1:size(part,1)
    patch('Faces',surfaces(sinpart==part(ii),:),'Vertices',nodes,'FaceColor',[38,255,255]/255,'LineWidth',1.0);
end
part=unique(einpart);
for ii=1:1:size(part,1)
    patch('Faces',elements((etype==223|etype==312)&einpart==part(ii),1:end-1),'Vertices',nodes,'FaceColor',[38,255,255]/255,'LineWidth',1.0);
end
%--------------------------------------------------------------------------
if size(nodes,2)==3
    view(-18,15)
end
axis equal
drawnow
end