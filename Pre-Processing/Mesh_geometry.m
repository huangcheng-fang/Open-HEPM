%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function New_Mesh=Mesh_geometry(geometry,meshSize,order)
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
model=createpde();model.Geometry=geometry;
generateMesh(model,"Hmin",meshSize(1,1),"Hmax",meshSize(1,2),'GeometricOrder',geomorder);
%==========================================================================
new_nodes=[];
new_Ninpart=[];
new_Nactivation=[];
new_elements=[];
new_Etype=[];
new_Einpart=[];
new_Eactivation=true(0,1);
pid=1;
temp_elements=model.Mesh.Elements'+size(new_nodes,1);
new_elements(end+1:end+size(temp_elements,1),1:size(temp_elements,2))=temp_elements;
new_nodes=[new_nodes;model.Mesh.Nodes'];
new_Ninpart=[new_Ninpart;repmat(pid,size(model.Mesh.Nodes,2),1)];
new_Nactivation=[new_Nactivation;true(size(model.Mesh.Nodes,2),1)];
new_Etype=[new_Etype;repmat(re_etype,size(temp_elements,1),1)];
new_Einpart=[new_Einpart;repmat(pid,size(temp_elements,1),1)];
new_Eactivation=[new_Eactivation;true(size(temp_elements,1),1)];
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
end