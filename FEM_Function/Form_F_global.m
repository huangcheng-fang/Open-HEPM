%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%%-------------------------------------------------------------------------
%Form global F vector
function F=Form_F_global(Load,Mesh,Set,DOF,max_Udof)
Key=Keywords();
nodes=Mesh.nodes;
elements=Mesh.elements;
Etype=Mesh.etype;
surfaces=Mesh.surfaces;
Stype=Mesh.stype;
Node_set=Set.node_set;
Surface_set=Set.surface_set;
Element_set=Set.element_set;
%==========================================================================
F=zeros(max_Udof,1);
for bi=1:numel(Load)
    setn=Load(bi).set(:);
    if isempty(setn);continue;end
    value=Load(bi).value(:);
    expression=Load(bi).expression;
    direction=Load(bi).direction;
    local_system=Load(bi).local_system;
    switch Load(bi).type
        case Key.concentrated_force
            set_obj=Merge_cell(Node_set,setn(:));
            F=Concentrated_force(F,nodes,set_obj,value,expression,direction,local_system,DOF);
        case Key.surface_traction
            set_obj=Merge_cell(Surface_set,setn(:));
            F=Surface_traction(F,nodes,surfaces,Stype,set_obj,value,expression,direction,local_system,DOF);
        case Key.pressure
            set_obj=Merge_cell(Surface_set,setn(:));
            F=Pressure(F,nodes,surfaces,Stype,set_obj,value,expression,DOF);
        case Key.body_force
            set_obj=Merge_cell(Element_set,setn(:));
            F=Body_force(F,nodes,elements,Etype,set_obj,value,expression,direction,local_system,DOF);
        case 'temp'
            set_obj=Merge_cell(Surface_set,setn(:));
            F=Pressure_temp(F,nodes,surfaces,Stype,set_obj,value,expression,DOF);
        otherwise
            error(['Unknown boundary condition type: ',Load(bi).type])
    end
end
%--------------------------------------------------------------------------
% time=toc;fprintf("External force calculation is completed:%fs\n",time);
end

function F=Concentrated_force(F,nodes,nid,value,expression,direction,local_system,DOF)
CF=value_or_expr(value,expression,nodes(nid,:));
CF=CF*local_system(direction,:);
dof=DOF(:,nid);
F(dof(:))=F(dof(:))+CF(:);
end

function F=Surface_traction(F,nodes,surfaces,Stype,set_obj,value,expression,direction,local_system,DOF)
% FI=zeros(maxdof,1);
[ip_tri,weight_tri]=Hammer_integral(2,3);
N_tri=Triangular_interpolation(ip_tri);
[ip_quad,weight_quad]=Gauss_Integral(2,3);
N_quad=Quadrilateral_interpolation(ip_quad);
dN_quad=Quadrilateral_interpolation_derivative(ip_quad);
[ip_tri_2order,weight_tri_2order]=Hammer_integral(2,3);
N_tri_2order=Triangular_quadratic_interpolation(ip_tri_2order);
dN_tri_2order=Triangular_quadratic_interpolation_derivative(ip_tri_2order);
for sk=1:1:size(set_obj,1)
    switch Stype(set_obj(sk))
        case 3
            surface_i=surfaces(set_obj(sk),1:3);
            F_dof=DOF(:,surface_i);
            node_coor=nodes(surface_i,:);
            detJ=norm(cross(node_coor(2,:)-node_coor(1,:),node_coor(3,:)-node_coor(1,:)));
            ip_global=N_tri*node_coor;
            f=value_or_expr(value,expression,ip_global);
            f=(f*local_system(direction,:)).';
            surface_F=f.*repmat((weight_tri.*detJ).',3,1)*N_tri;
            F(F_dof)=F(F_dof)+surface_F;
        case 4
            surface_i=surfaces(set_obj(sk),1:4);
            F_dof=DOF(:,surface_i);
            node_coor=nodes(surface_i,:);
            Normal_vector=Get_normal_vector(node_coor);
            R=vector_to_system(Normal_vector);
            node_coorp=node_coor*R(1:2,:).';
            detJ=zeros(size(dN_quad,3),1);
            for j=1:1:size(dN_quad,3)
                detJ(j)=abs(det(dN_quad(:,:,j)*node_coorp));
            end
            ip_global=N_quad*node_coor;
            f=value_or_expr(value,expression,ip_global);
            f=(f*local_system(direction,:)).';
            surface_F=f.*repmat((weight_quad.*detJ).',3,1)*N_quad;
            F(F_dof)=F(F_dof)+surface_F;
        case 6
            surface_i=surfaces(set_obj(sk),1:6);
            F_dof=DOF(:,surface_i);
            node_coor=nodes(surface_i,:);
            Normal_vector=Get_normal_vector(node_coor([1,3,5],:));
            R=vector_to_system(Normal_vector);
            node_coorp=node_coor*R(1:2,:).';
            detJ=zeros(size(dN_tri_2order,3),1);
            for j=1:1:size(dN_tri_2order,3)
                detJ(j)=abs(det(dN_tri_2order(:,:,j)*node_coorp));
            end
            ip_global=N_tri_2order*node_coor;
            f=value_or_expr(value,expression,ip_global);
            f=(f*local_system(direction,:)).';
            surface_F=f.*repmat((weight_tri_2order.*detJ).',3,1)*N_tri_2order;
            F(F_dof)=F(F_dof)+surface_F;
    end
end
end

function F=Pressure(F,nodes,surfaces,Stype,set_obj,value,expression,DOF)
% FI=zeros(maxdof,1);
[ip_tri,weight_tri]=Hammer_integral(2,3);
SN_tri=Triangular_interpolation(ip_tri);
[ip_quad,weight_quad]=Gauss_Integral(2,3);
N_quad=Quadrilateral_interpolation(ip_quad);
dN_quad=Quadrilateral_interpolation_derivative(ip_quad);
for sk=1:1:size(set_obj,1)
    switch Stype(set_obj(sk))
        case 3
            surface_i=surfaces(set_obj(sk),1:3);
            F_dof=DOF(:,surface_i);
            node_coor=nodes(surface_i,:);

            normal_vector=-cross(node_coor(2,:)-node_coor(1,:),node_coor(3,:)-node_coor(1,:));
            detJ=norm(normal_vector);normal_vector=normal_vector./detJ;

            ip_global=SN_tri*node_coor;
            f=value_or_expr(value,expression,ip_global);
            f=(f*normal_vector).';
            surface_F=f.*repmat((weight_tri.*detJ).',3,1)*SN_tri;
            F(F_dof)=F(F_dof)+surface_F;
        case 4
            surface_i=surfaces(set_obj(sk),1:4);
            F_dof=DOF(:,surface_i);
            node_coor=nodes(surface_i,:);
            Normal_vector=-Get_normal_vector(node_coor);
            R=vector_to_system(Normal_vector);
            node_coorp=node_coor*R(1:2,:).';
            detJ=zeros(size(dN_quad,3),1);
            for j=1:1:size(dN_quad,3)
                detJ(j)=abs(det(dN_quad(:,:,j)*node_coorp));
            end
            ip_global=N_quad*node_coor;
            f=value_or_expr(value,expression,ip_global);
            f=(f*Normal_vector).';
            surface_F=f.*repmat((weight_quad.*detJ).',3,1)*N_quad;
            F(F_dof)=F(F_dof)+surface_F;
        otherwise
            error('Unsurpported surface type')
    end
end
end

function F=Pressure_temp(F,nodes,surfaces,Stype,set_obj,value,expression,DOF)
% FI=zeros(maxdof,1);
[ip_tri,weight_tri]=Hammer_integral(2,3);
SN_tri=Triangular_interpolation(ip_tri);
ref_point=nodes(value(2),:);
for sk=1:1:size(set_obj,1)
    switch Stype(set_obj(sk))
        case 3
            surface_i=surfaces(set_obj(sk),1:3);
            F_dof=DOF(:,surface_i);
            node_coor=nodes(surface_i,:);

            normal_vector=-cross(node_coor(2,:)-node_coor(1,:),node_coor(3,:)-node_coor(1,:));
            detJ=norm(normal_vector);normal_vector=normal_vector./detJ;

            ip_global=SN_tri*node_coor;

            v=cross(normal_vector,[0,1,0]);
            f=(ip_global-ref_point)*v.'*value(1);

%             f=value_or_expr(value,expression,ip_global);
            f=(f*normal_vector).';
            surface_F=f.*repmat((weight_tri.*detJ).',3,1)*SN_tri;
            F(F_dof)=F(F_dof)+surface_F;
        otherwise
            error('Unsurpported surface type')
    end
end
end

function F=Body_force(F,nodes,elements,Etype,set_obj,value,expression,direction,local_system,DOF)
% FI=zeros(maxdof,1);
[ip_tet,weight_tet]=Hammer_integral(3,2);
N_tet=Tetrahedron_interpolation(ip_tet);
[ip_hex,weight_hex]=Gauss_Integral(3,3);
N_hex=Hexahedron_interpolation(ip_hex);
dN_hex=Hexahedron_interpolation_derivative(ip_hex);
for ek=1:1:size(set_obj,1)
    switch Etype(set_obj(ek))
        case 334
            element=elements(set_obj(ek),1:4);
            F_dof=DOF(:,element);
            node_coor=nodes(element,1:3);
            Jac=node_coor(2:end,:)-repmat(node_coor(1,:),3,1);
            detJ=det(Jac);
            ip_global=N_tet*node_coor;
            f=value_or_expr(value,expression,ip_global);
            f=(f*local_system(direction,:)).';
            Body_F=f.*repmat((weight_tet.*detJ).',3,1)*N_tet;
            F(F_dof)=F(F_dof)+Body_F;
        case 338
            element=elements(set_obj(ek),1:8);
            F_dof=DOF(:,element);
            node_coor=nodes(element,1:3);
            detJ=zeros(size(dN_hex,3),1);
            for j=1:1:size(dN_hex,3)
                detJ(j)=abs(det(dN_hex(:,:,j)*node_coor));
            end
            ip_global=N_hex*node_coor;
            f=value_or_expr(value,expression,ip_global);
            f=(f*local_system(direction,:)).';
            Body_F=f.*repmat((weight_hex.*detJ).',3,1)*N_hex;
            F(F_dof)=F(F_dof)+Body_F;
        otherwise
            error('Unsurpported element type')
    end
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