%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%%-------------------------------------------------------------------------
%Form global F vector
function F=Form_fluid_flux_global(Flux_condition,Mesh_Topology,Set)
Key=Keywords();
nodes=Mesh_Topology.nodes;
elements=Mesh_Topology.elements;
surfaces=Mesh_Topology.surfaces;
Node_set=Set.node_set;
Surface_set=Set.surface_set;
Element_set=Set.element_set;
%==========================================================================
F=zeros(size(nodes,1),1);
for bi=1:numel(Flux_condition)
    setn=Flux_condition(bi).set(:);
    if isempty(setn);continue;end
    value=Flux_condition(bi).value(:);
    expression=Flux_condition(bi).expression;
    switch Flux_condition(bi).type
        case Key.nodal_flux
            set_obj=Merge_cell(Node_set,setn(:));
            F=Nodal_flux(F,nodes,set_obj,value,expression);
        case Key.surface_flux
            set_obj=Merge_cell(Surface_set,setn(:));
            F=Surface_flux(F,nodes,surfaces,set_obj,value,expression);
        case Key.body_flux
            set_obj=Merge_cell(Element_set,setn(:));
            F=Body_flux(F,nodes,elements,set_obj,value,expression);
        otherwise
            error(['Unknown boundary condition type: ',Flux_condition(bi).type])
    end
end
%--------------------------------------------------------------------------
% time=toc;fprintf("External force calculation is completed:%fs\n",time);
end

function F=Nodal_flux(F,nodes,nid,value,expression)
% FI=zeros(maxdof,1);
CF=value_or_expr(value,expression,nodes(nid,:));
dof=nid;
F(dof)=F(dof)+CF;
end

function F=Surface_flux(F,nodes,surfaces,set_obj,value,expression)
% FI=zeros(maxdof,1);
[ip,weight]=Triangular_integral(2);
SN=Triangular_interpolation(ip);
for sk=1:1:size(set_obj,1)
    surface_i=surfaces(set_obj(sk),1:3);
    F_dof=surface_i;
    node_coor=nodes(surface_i,:);
    
    normal_area=cross(node_coor(2,:)-node_coor(1,:),node_coor(3,:)-node_coor(1,:));
    detJ=norm(normal_area);
    
    ip_global=SN*node_coor;
    f=value_or_expr(value,expression,ip_global);
    surface_F=(f.*weight.*detJ).'*SN;
    F(F_dof)=F(F_dof)+surface_F;
end
end

function F=Body_flux(F,nodes,elements,set_obj,value,expression)
% FI=zeros(maxdof,1);
[ip,weight]=Tetrahedron_integral(2);
SN=Tetrahedron_interpolation(ip);
for ek=1:1:size(set_obj,1)
    element=elements(set_obj(ek),1:4);
    F_dof=element;
    node_coor=nodes(element,1:3);
    
    Jac=node_coor(2:end,:)-repmat(node_coor(1,:),3,1);
    detJ=det(Jac);
    
    ip_global=SN*node_coor;
    f=value_or_expr(value,expression,ip_global);
    Body_F=(f.*weight.*detJ).'*SN;
    F(F_dof)=F(F_dof)+Body_F;
end
end

function V=value_or_expr(value,expression,nodes)
if numel(expression)~=numel(value)
    error('The size of Expression should be same with Value');
end
x=nodes(:,1);y=nodes(:,2);
V=zeros(numel(x),numel(expression));
if size(nodes,2)==3;z=nodes(:,3);else;z=zeros(numel(x),1);end
for epi=1:1:numel(expression)
    if isempty(expression{epi})
        V(:,epi)=repmat(value(epi),size(nodes,1),1);
    else
        V(:,epi)=real(Simple_calculator(expression{epi},x,y,z));
    end
end
end