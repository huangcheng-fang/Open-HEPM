%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%By Fang Huangcheng @PolyU
%Implicit return algorithm; Backward Euler Method
%Last update @2024/1/31
%Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
%Researchgate:www.researchgate.net/profile/Huangcheng-Fang
function New_Mesh=Copy_Remesh_domain_new(Mesh,rangepart)
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
Stype=Mesh.stype;
%==========================================================================
[surfaces,Sinpart]=split_quadrilateral_surface(surfaces,Sinpart,Stype);
Nconnection=Get_node_connection(elements);
Nactivation=deactive_too_close_points(Nconnection,Nactivation,nodes,elements,surfaces,0.1);
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
    flag=1;
    while true
    Index=zeros(size(nodes,1),1);
    node_index=find(Ninpart==partid(ii)&Nactivation);
    part_nodes=nodes(node_index,:);
    Index(node_index)=1:1:numel(node_index);
    part_surface=Index(surfaces(Sinpart==partid(ii),:));
    %--------------------------------------------------------------------------
%     DT=delaunayTriangulation(part_nodes);
    DT=constrainedDelaunayTetGen(part_nodes,part_surface);
    part_elements=DT.ConnectivityList;
    %----------------------------remove outer element--------------------------
    RMloc=find_outer_element(part_nodes,part_elements,part_surface);
    part_elements(RMloc,:)=[];
    isdistored=find_distored_element(part_elements,part_nodes,0.1);
    [part_elements,true_distored_nid]=remove_internal_slivers_element(part_nodes,part_elements,find(isdistored));
    true_distored_nid(ismember(true_distored_nid,part_surface))=[];
    Nactivation(node_index(true_distored_nid))=false;
    flag=flag+1;
    if isempty(true_distored_nid);break;end
    if flag>100
        error('in remesh')
    end
    end
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

function Nactivation=deactive_too_close_points(Nconnection,Nactivation,nodes,elements,surfaces,h_tol)
issurface=false(size(Nconnection,1),1);
nid=surfaces;nid(nid==0)=[];
issurface(nid)=true;
for ni=1:1:size(Nconnection,1)
if Nactivation(ni)&&~issurface(ni)
    eid=Nconnection(ni,1:Nconnection(ni,end));
    nid=unique(elements(eid,:));
    nid(nid==ni)=[];
    loc=vecnorm(nodes(nid,:)-nodes(ni,:),2,2)<h_tol;
    Nactivation(nid(loc))=false;
end
end
Nactivation(issurface)=true;
end


function RMloc=find_outer_element(part_nodes,part_elements,part_surface)
boundary_element_id=find(all(ismember(part_elements,part_surface),2));
boundary_elements=part_elements(boundary_element_id,:);
isdistored=find_distored_element(boundary_elements,part_nodes,1e-1);
RMloc=boundary_element_id(isdistored,:);
boundary_element_id(isdistored,:)=[];
boundary_elements(isdistored,:)=[];
boundary_first_node_id=unique(boundary_elements(:,1));
%=======================find active surface================================
surface_connection=Get_node_connection(part_surface);
active_surface=unique(surface_connection(boundary_first_node_id,1:end-1));
active_surface(active_surface==0)=[];
part_surface=part_surface(active_surface,:);
surface_connection=Get_node_connection(part_surface);
%==========================surface normal==================================
line1=part_nodes(part_surface(:,2),:)-part_nodes(part_surface(:,1),:);
line2=part_nodes(part_surface(:,3),:)-part_nodes(part_surface(:,1),:);
surface_normal=cross(line1,line2,2);
surface_area_sqrt=sqrt(sum(surface_normal.^2,2));
surface_normal=surface_area_sqrt.\surface_normal;
surface_area_sqrt=sqrt(surface_area_sqrt);
%====================remove outer element==================================
element_center=(part_nodes(boundary_elements(:,1),:)+part_nodes(boundary_elements(:,2),:)+part_nodes(boundary_elements(:,3),:)+part_nodes(boundary_elements(:,4),:))/4;
remove_flag=false(size(boundary_elements,1),1);
for ei=1:1:size(boundary_elements,1)

    nid=boundary_elements(ei,1);
    sid=unique(surface_connection(nid,1:surface_connection(nid,end)));
    snid=part_surface(sid(1),:);
    ref_point=surface_normal(sid(1),:)*(1e-8*surface_area_sqrt(sid(1)))+mean(part_nodes(snid,:),1)-part_nodes(nid,:);
    emid_point=element_center(ei,:)-part_nodes(nid,:);
    flag=0;
    for i=1:1:numel(sid)
        % hold on;patch('Faces',[1,2,3],'Vertices',part_nodes(part_surface(sid(i),:),1:3),'FaceColor',[rand(1),rand(1),rand(1)],'LineWidth',0.01);
        d1=(surface_normal(sid(i),:)*ref_point.');
        d2=(surface_normal(sid(i),:)*emid_point.');

        if abs(d2)/surface_area_sqrt(sid(1))<1e-12
            flag=0;break
        end
        if d1*d2>=0
            continue
        end
        d1=abs(d1);d2=abs(d2);

        intersection=(ref_point*d2+emid_point*d1)/(d1+d2);
        snid=part_surface(sid(i),:);
        loc=snid~=nid;
        surface_line=part_nodes(snid(loc),:)-part_nodes(nid,:);

        surface_line=sqrt(sum(surface_line.^2,2)).\surface_line;
        intersection=norm(intersection).\intersection;
        if all(surface_line(1,:)*surface_line(2,:).'<surface_line*intersection')
            flag=flag+1;
        end
    end
    if mod(flag,2)==0
        remove_flag(ei)=true;
    else
        % hold on;plot_element_temp(ei,boundary_elements,part_nodes,repmat(334,100000,1))
    end
end
RMloc=[RMloc;boundary_element_id(remove_flag)];
end

function isdistored=find_distored_element(tet_elements,nodes,tol)
isdistored=true(size(tet_elements,1),1);
for ei=1:1:size(tet_elements,1)
    nid=tet_elements(ei,1:4);
    element_nodes=nodes(nid,1:3);
    Jac=element_nodes(2:end,:)-repmat(element_nodes(1,:),3,1);
    % isdistored(ei)=det(Jac)/sum(Jac.^2,'all')<tol;
    a(1)=norm(cross(Jac(1,:),Jac(2,:)));
    a(2)=norm(cross(Jac(2,:),Jac(3,:)));
    a(3)=norm(cross(Jac(3,:),Jac(1,:)));
    h=det(Jac)/max(a);
    l=max(vecnorm(Jac'));
    isdistored(ei)=h/l<tol;
end
end

function part_elements=remove_internal_thin_element(part_nodes,part_elements)
isdistored=find_distored_element(part_elements,part_nodes,1e-2);
dist_element_list=find(isdistored);
nconnection=Get_node_connection(part_elements);

for i=1:1:numel(dist_element_list)
dist_eid=dist_element_list(i);
dist_enid=part_elements(dist_eid,:);
connected_eid=unique(nconnection(dist_enid,1:end-1));
connected_eid(connected_eid==0)=[];
loc=sum(ismember(part_elements(connected_eid,:),dist_enid),2)>2;
connected_eid=connected_eid(loc);
connected_eid(connected_eid==dist_eid)=[];
if numel(connected_eid)~=4
    figure(100);clf(100)
    tetramesh(part_elements([connected_eid;dist_eid],:),part_nodes)
end
continue
if numel(connected_eid)~=4;error('This is a reserved check for potential bugs');end

connected_elements=part_elements(connected_eid,:);
loc=ismember(connected_elements,dist_enid);
nid=unique(connected_elements(~loc));
loc=any(ismember(connected_elements,nid(1)),2);
temp=connected_elements(loc,[2,1,3,4]);
temp(temp==nid(1))=nid(2);
connected_elements(~loc,:)=temp;
part_elements(connected_eid,:)=connected_elements;
end

part_elements(isdistored,:)=[];
end


function [surfaces,Sinpart]=split_quadrilateral_surface(surfaces,Sinpart,Stype)
for stypei=4:1:max(Stype,[],'all')
    loc=Stype==stypei;

    if any(loc,1)
        temp_Sinpart=Sinpart(loc,:);
        temp_Stype=Stype(loc,:);
        temp_surfaces=surfaces(loc,:);

        Sinpart(loc,:)=[];
        Stype(loc,:)=[];
        surfaces(loc,:)=[];

        for i=1:1:stypei-2
            Sinpart=[Sinpart;temp_Sinpart];
            Stype=[Stype;temp_Stype];
            surfaces(end+1:end+size(temp_surfaces,1),1:3)=temp_surfaces(:,[1,i+1,i+2]);
        end
    end
end

surfaces=surfaces(:,1:3);
end