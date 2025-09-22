%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Plot_set(Model,n,type,setid)
nodes=Model.Mesh.nodes;
elements=Model.Mesh.elements;
surfaces=Model.Mesh.surfaces;
etype=Model.Mesh.etype;
Set=Model.Set;
%---------------------------identify type----------------------------------
nset=[];eset=[];sset=[];
switch type
    case 'node'
        nset=setid;
    case 'element'
        eset=setid;
    case 'surface'
        sset=setid;
    otherwise
        error('Unkonwn type')
end
%------------------------------plot node-----------------------------------
figure(n);hold on;
for i=1:1:numel(nset)
    NID=Set.node_set{nset(i)}(:);
    if size(nodes,2)==2
        plot(nodes(NID,1),nodes(NID,2),'.','MarkerSize',12)
    else
        plot3(nodes(NID,1),nodes(NID,2),nodes(NID,3),'.','MarkerSize',20);
    end
end
%------------------------------plot surface--------------------------------
for i=1:1:numel(sset)
    SID=Set.surface_set{sset(i)}(:);
    if size(nodes,2)==2
        patch('Faces',surfaces(SID,:),'Vertices',nodes,'FaceColor',[rand(1),rand(1),rand(1)],'LineWidth',1,'EdgeColor','red');
    else
        patch('Faces',surfaces(SID,:),'Vertices',nodes,'FaceColor',[rand(1),rand(1),rand(1)],'LineWidth',0.001);
    end
end
%----------------------------plot element----------------------------------
for i=1:1:numel(eset)
    EID=Set.element_set{eset(i)}(:);
    plot_element(EID(:),elements,nodes,etype)
end
%-------------------------------3D view------------------------------------
axis equal;
% if size(nodes,2)==3
%     view(-18,15)
% end
drawnow;
end


function plot_element(EID,elements,nodes,type)
faces=nan(size(EID,1)*6,4);
for ii=1:1:size(EID,1)
    switch type(EID(ii))
        case 223
            E=elements(EID(ii),1:3);
            faces(ii*6,1:3)=E;
        case 224
            E=elements(EID(ii),1:4);
            faces(ii*6,1:4)=E;
        case 334
            face=[2,1,3;1,2,4;2,3,4;3,1,4];
            E=elements(EID(ii),1:4);
            faces(ii*6-3:ii*6,1:3)=E(face);
        case {338,3381}
            face=[4,3,2,1;5,6,7,8;6,5,1,2;2,3,7,6;8,7,3,4;1,5,8,4];
            E=elements(EID(ii),1:8);
            faces(ii*6-5:ii*6,1:4)=E(face);
        case 323
            face=[1,2,3];
            E=elements(EID(ii),1:3);
            faces(ii*6,1:3)=E(face);
        case 324
            face=[1,2,3,4];
            E=elements(EID(ii),1:4);
            faces(ii*6,1:4)=E(face);
        case 324.4
            face=[1,2,3,4];
            E=elements(EID(ii),1:4);
            faces(ii*6,1:4)=E(face);
        case 312
            E=elements(EID(ii),1:2);
            faces(ii*6,1:2)=E;
        case 312.02
            E=elements(EID(ii),1:2);
            faces(ii*6,1:2)=E;
        case 312.2
            E=elements(EID(ii),1:2);
            faces(ii*6,1:2)=E;
        otherwise
            error('Unknown element type');
    end
end
patch('Faces',faces,'Vertices',nodes,'FaceColor',[rand(1),rand(1),rand(1)],'LineWidth',0.01);
end