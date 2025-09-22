%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function New_Mesh=Remesh_geometry(Mesh,rangePart,meshSize,faceSize,order,addNodes)
Mesh=Get_surface(Mesh);
nodes=Mesh.nodes;
Ninpart=Mesh.ninpart;
Nactivation=Mesh.nactivation;
elements=Mesh.elements;
Etype=Mesh.etype;
Einpart=Mesh.einpart;
Eactivation=Mesh.eactivation;
surfaces=Mesh.surfaces;
Sinpart=Mesh.sinpart;
Stype=Mesh.stype;
particles=Mesh.particles;
pinpart=Mesh.pinpart;
%==========================================================================
if order==1
    geomorder='linear';
    re_etype=334;
elseif order==2
    geomorder='quadratic';
    re_etype=3310;
else
    error('unsurpported GeometricOrder')
end
%==========================================================================
if size(surfaces,2)<10;surfaces(:,10)=0;end
rangePart=unique(rangePart);
for i=1:1:numel(rangePart)
    loc=Sinpart==rangePart(i);
    part_surfaces=surfaces(loc,:);
    part_stype=Stype(loc);
    loc3=part_stype==3;
    loc4=part_stype==4;
    part_surfaces=[part_surfaces(loc3,1:3);
                   part_surfaces(loc4,1:3);
                   part_surfaces(loc4,[3,4,1])];
    model(i)=createpde();
    geometryFromMesh(model(i),nodes',part_surfaces');
    if ~isempty(addNodes{i})
        addVertex(model(i).Geometry,"Coordinates",addNodes{i});
    end
    for ii=1:2:numel(faceSize)
        faceSize{ii}=nearestFace(model(i).Geometry,faceSize{ii});
    end
    generateMesh(model(i),"Hmin",meshSize(i,1),"Hmax",meshSize(i,2),"Hface",faceSize,"Hgrad",1.1,'GeometricOrder',geomorder);
end
%==========================================================================
elements(elements==0)=nan;
new_nodes=[];
new_Ninpart=[];
new_Nactivation=[];
new_elements=[];
new_Etype=[];
new_Einpart=[];
new_Eactivation=true(0,1);
for mi=1:1:numel(rangePart)
    pid=rangePart(mi);
    temp_elements=model(mi).Mesh.Elements'+size(new_nodes,1);
    new_elements(end+1:end+size(temp_elements,1),1:size(temp_elements,2))=temp_elements;
    new_nodes=[new_nodes;model(mi).Mesh.Nodes'];
    new_Ninpart=[new_Ninpart;repmat(pid,size(model(mi).Mesh.Nodes,2),1)];
    new_Nactivation=[new_Nactivation;true(size(model(mi).Mesh.Nodes,2),1)];
    new_Etype=[new_Etype;repmat(re_etype,size(temp_elements,1),1)];
    new_Einpart=[new_Einpart;repmat(pid,size(temp_elements,1),1)];
    new_Eactivation=[new_Eactivation;true(size(temp_elements,1),1)];
end
%==========================================================================
eid=find(~ismember(Einpart,rangePart));
temp_elements=elements(eid,:);
nid=unique(temp_elements(:));
nid(isnan(nid))=[];
Index=zeros(max(nid),1);
Index(nid,1)=size(new_nodes,1)+1:size(new_nodes,1)+numel(nid);
Index=[Index;nan];
temp_elements(isnan(temp_elements))=size(Index,1);
new_elements(end+1:end+size(temp_elements,1),1:size(temp_elements,2))=Index(temp_elements);
new_nodes=[new_nodes;nodes(nid,:)];
new_Ninpart=[new_Ninpart;Ninpart(nid)];
new_Nactivation=[new_Nactivation;Nactivation(nid)];
new_Etype=[new_Etype;Etype(eid)];
new_Einpart=[new_Einpart;Einpart(eid)];
new_Eactivation=[new_Eactivation;Eactivation(eid)];
%==========================================================================
new_elements(new_elements==0)=nan;
New_Mesh=Creat_new_mesh();
New_Mesh.nodes=new_nodes;
New_Mesh.positions=new_nodes;
New_Mesh.ninpart=new_Ninpart;
New_Mesh.nactivation=new_Nactivation;
New_Mesh.elements=new_elements;
New_Mesh.etype=new_Etype;
New_Mesh.einpart=new_Einpart;
New_Mesh.ecenter=Get_geometry_center(new_elements,new_nodes);
New_Mesh.eactivation=new_Eactivation;
New_Mesh=Get_surface(New_Mesh);
New_Mesh.particles=particles;
New_Mesh.pinpart=pinpart;
end