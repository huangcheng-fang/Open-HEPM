%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%By Fang Huangcheng @PolyU
%Email: huangcheng.fang@polyu.edu.hk; valy_f@bjtu.edu.cn
function New_Mesh=Remesh_domain(Mesh,rangepart)
tic;
nodes=Mesh.nodes;
positions=Mesh.positions;
Ninpart=Mesh.ninpart;
Nactivation=Mesh.nactivation;
elements=Mesh.elements;
Einpart=Mesh.einpart;
Etype=Mesh.etype;
Eactivation=Mesh.eactivation;
surfaces=Mesh.surfaces;
Sinpart=Mesh.sinpart;
surfaces(isnan(surfaces))=size(nodes,1)+1;
%==========================main function===================================
partid=unique(Ninpart);
container=Create_parallel_container(numel(partid));
%--------------------------------------------------------------------------
for ii=1:numel(partid)
if ~ismember(partid(ii),rangepart)
    loc=Einpart==partid(ii);
    container(ii).c1=elements(loc,:);
    container(ii).c2=repmat(partid(ii),size(container(ii).c1,1),1);
    container(ii).c3=Etype(loc);
    container(ii).c4=Eactivation(loc);
    continue;
end
Index=zeros(size(nodes,1)+1,1);Index(end)=nan;
node_index=find(Ninpart==partid(ii)&Nactivation);
part_nodes=nodes(node_index,:);
Index(node_index)=1:1:numel(node_index);
part_surface=Index(surfaces(Sinpart==partid(ii),:));
part_elements0=elements(Einpart==partid(ii)&Eactivation,:);
part_elements0(:,isnan(part_elements0(1,:)))=[];
part_elements0=Index(part_elements0);
%--------------------------------------------------------------------------
DT=delaunayTriangulation(part_nodes);
part_elements=DT.ConnectivityList;
%----------------------------remove outer element--------------------------
RMloc=remove_outer_element(part_nodes,part_elements,part_elements0,part_surface);
part_elements(RMloc,:)=[];
%--------------------------------------------------------------------------
container(ii).c1=node_index(part_elements);
container(ii).c2=repmat(partid(ii),size(part_elements,1),1);
container(ii).c3=repmat(334,size(part_elements,1),1);
container(ii).c4=true(size(part_elements,1),1);
end
[elements,Einpart,Etype,Eactivation]=Merge_container(container,1);elements(elements==0)=nan;
%==============================Output======================================
New_Mesh=Creat_new_mesh();
New_Mesh.nodes=nodes;
New_Mesh.positions=positions;
New_Mesh.ninpart=Ninpart;
New_Mesh.nactivation=Nactivation;
% New_Mesh.nconnection=Get_node_connection(size(nodes,1),elements,Etype);
New_Mesh.elements=elements;
New_Mesh.etype=Etype;
New_Mesh.einpart=Einpart;
New_Mesh.eactivation=logical(Eactivation);
New_Mesh.surfaces=Mesh.surfaces;
New_Mesh.sinfacet=Mesh.sinfacet;
New_Mesh=Get_surface(New_Mesh);
New_Mesh=Update_mesh_geometry(New_Mesh);
%==============================disp========================================
time=toc;fprintf("Computational domain has been remeshed: %fs\n",time);
end


function RMloc=remove_outer_element(part_nodes,part_elements,part_elements0,part_surface)
RMloc=false(size(part_elements,1),1);
removeloc1=find(all(ismember(part_elements,part_surface),2));
po_elements=part_elements(removeloc1,:);
%--------------------------------------------------------------------------
if size(part_elements0,2)==8
part_elements0=[part_elements0(:,[7 6 8 3]);
part_elements0(:,[5 8 6 3]);
part_elements0(:,[5 6 2 3]);
part_elements0(:,[5 4 8 3]);
part_elements0(:,[5 2 4 3]);
part_elements0(:,[5 2 1 4])];
end
nconnection=Get_node_connection(part_elements0);
%--------------------------------------------------------------------------
removeloc2=remove_distored_element(po_elements,part_nodes);
for ei=1:1:size(po_elements,1)
    if removeloc2(ei);continue;end
    emid=mean(part_nodes(po_elements(ei,:),:),1);
    eid=unique(nconnection(po_elements(ei,:),:));
    eid(eid==0)=[];
    flag=true;
    for ii=1:1:numel(eid)
        enid=part_elements0(eid(ii),1:4);
        ip=Tetrahedron_inverse_interpolation(part_nodes(enid,1:3),emid);
        if min(ip)>=-1e-8&&max(ip)<=1.000000001&&sum(ip)<=1.00000001
            flag=false;break
        end
    end
    removeloc2(ei)=flag;
end
%--------------------------------------------------------------------------
RMloc(removeloc1(removeloc2))=true;
end


function isdistored=remove_distored_element(tet_elements,nodes)
isdistored=true(size(tet_elements,1),1);
for ei=1:1:size(tet_elements,1)
    nid=tet_elements(ei,1:4);
    element_nodes=nodes(nid,1:3);
    Jac=element_nodes(2:end,:)-repmat(element_nodes(1,:),3,1);
    isdistored(ei)=det(Jac)/sum(Jac.^2,'all')<1e-4;
end
end
