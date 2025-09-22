%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FE=FEM_Element(Mesh,MID)
nodes=Mesh.nodes;
elements=Mesh.elements;
Etype=Mesh.etype;
%===================4-node tetrahedron element===============
tet4loc=find(Etype==334);
if ~isempty(tet4loc)
    tet4e=elements(tet4loc,1:4);
    Tet4Element=Tetrahedron_element(tet4e,nodes);
    Tet4Element.eid=tet4loc;
    [ip,~]=Tetrahedron_integral(1);
    Tet4Element.N_matrix=Tetrahedron_interpolation(ip);
    FE(1,1)=Tet4Element;
end
%===========================8-node Hexahedron element======================
hex8loc=find(Etype==338);
if ~isempty(hex8loc)
    IntegralOder=2;
    hex8e=elements(hex8loc,1:8);
    Hex8Element=Hexahedron_element(hex8e,nodes,IntegralOder);
    Hex8Element.eid=hex8loc;
    [ip,~]=Gauss_Integral(3,IntegralOder);
    Hex8Element.N_matrix=permute(Hexahedron_interpolation(ip),[3,2,1]);
    FE(2,1)=Hex8Element;
end
%===================3-node triangle element=================
tri3loc=find(Etype==323);
if ~isempty(tri3loc)
    tri3e=elements(tri3loc,1:3);
    Tri3Element=Triangular_Element_323(tri3e,nodes);
    Tri3Element.eid=tri3loc;
    [ip,~]=Triangular_integral(1);
    Tri3Element.N_matrix=permute(Triangular_interpolation(ip),[3,2,1]);
    FE(3,1)=Tri3Element;
end
%===================4-node Quadrilateral element=================
quad324loc=find(Etype==324);
if ~isempty(quad324loc)
    IntegralOder=2;
    quad324e=elements(quad324loc,1:4);
    Quad324Element=Quadrilateral_Element_324(quad324e,nodes,IntegralOder);
    Quad324Element.eid=quad324loc;
    ip=Gauss_Integral(2,IntegralOder);
    Quad324Element.N_matrix=permute(Quadrilateral_interpolation(ip),[3,2,1]);
    FE(4,1)=Quad324Element;
end
%===================2-node link element=================
link3loc=find(Etype==312);
if ~isempty(link3loc)
    IntegralOder=1;
    link3e=elements(link3loc,1:2);
    Link3Element=Link_Element_312(link3e,nodes,IntegralOder);
    Link3Element.eid=link3loc;
    ip=Gauss_Integral(1,IntegralOder);
    Link3Element.N_matrix=permute(Line_interpolation(ip),[3,2,1]);
    FE(5,1)=Link3Element;
end
%===================8-node incompressible Hexahedron element===============
hex8loc=find(Etype==3381);
if ~isempty(hex8loc)
    IntegralOder=2;
    hex8e=elements(hex8loc,1:8);
    Hex8Element=Hexahedron_element_incompressible_3381(hex8e,nodes,IntegralOder);
    Hex8Element.eid=hex8loc;
    [ip,~]=Gauss_Integral(3,IntegralOder);
    Hex8Element.N_matrix=permute(Hexahedron_interpolation(ip),[3,2,1]);
    FE(5,1)=Hex8Element;
end
%=========================Assign_material====================
FE=Assign_FE_material(FE,MID);
end

function FEM=Assign_FE_material(FE,MID)
FEM(numel(FE)*max(MID),1).type=zeros(1,1);
ii=1;
for fi=1:1:numel(FE)
    mis=MID(FE(fi).eid);
    mino=unique(mis);
    if isempty(mis);continue;end
    for mi=1:1:numel(mino)
        loc=mis==mino(mi);
        FEM(ii).type=FE(fi).type;
        FEM(ii).R_matrix=FE(fi).R_matrix(:,:,:,loc);
        FEM(ii).B_matrix=FE(fi).B_matrix(:,:,:,loc);
        FEM(ii).FB_matrix=FE(fi).FB_matrix(:,:,:,loc);
        FEM(ii).detJac=FE(fi).detJac(:,:,:,loc);
        FEM(ii).dof=FE(fi).dof(:,loc);
        FEM(ii).eid=FE(fi).eid(loc,1);
        FEM(ii).N_matrix=FE(fi).N_matrix;
        FEM(ii).materialID=mino(mi);
        ii=ii+1;
    end
end
FEM(ii:end)=[];
end

% [Econstant,map,usable_type]=Element_constant_list();
% for ti=1:1:numel(usable_type)
%     loc=find(Etype==usable_type(ti));
%     if ~isempty(loc)
%         Enid=elements(loc,1:Econstant(ti).node_num);
%         IntegralOder=Econstant(ti).int_order;
%         switch usable_type(ti)
%             case 334
%                 Element=Tetrahedron_element(Enid,nodes);
%             case 338
%                 Element=Hexahedron_element(Enid,nodes,IntegralOder);
%             case 323
%                 Element=Triangular_Element_323(Enid,nodes);
%             case 324
%                 Element=Quadrilateral_Element_324(Enid,nodes,IntegralOder);
%             case 312
%                 Element=Link_Element_312(Enid,nodes,IntegralOder);
%             case 3381    
%                 Element=Hexahedron_element_incompressible_3381(Enid,nodes,IntegralOder);
%             otherwise
%                 error('Unknown element type')
%         end
%         Element.eid=loc;
%         Element.N_matrix=permute(Econstant(ti).N_matrix,[3,2,1]);
%         FE(ti,1)=Element;
%     end
% end