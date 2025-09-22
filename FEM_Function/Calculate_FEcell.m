%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FEcell=Calculate_FEcell(FEcell,Mesh,field_dof_index)
tic;
nodes=Mesh.nodes;
elements=Mesh.elements;
Eactivation=Mesh.eactivation;
[Econstant,Etype_list]=Element_constant_list();
%==============================main function===============================
for fi=1:1:numel(FEcell.Element_id)
    type=FEcell.Type{fi};
    loc=FEcell.Element_id{fi};
    loc=loc(Eactivation(loc));
    ti=find(Etype_list==type,1);
    if ~isempty(loc)
        Enid=elements(loc,1:Econstant(ti).node_num);
        IntegralOrder=Econstant(ti).int_order;
        switch type
            case 334
                Element=Tetrahedron_element(Enid,nodes,IntegralOrder);
            case 338
                Element=Hexahedron_element(Enid,nodes,IntegralOrder);
            case 323
                Element=Triangular_Element_323(Enid,nodes);
            case 324
                Element=Quadrilateral_Element_324(Enid,nodes,IntegralOrder);
            case 312
                Element=Link_Element_312(Enid,nodes,IntegralOrder);
            case 3381    
                Element=Hexahedron_element_incompressible_3381(Enid,nodes,IntegralOrder);
            case 3310
                Element=Tetrahedron_quadratic_element(Enid,nodes,IntegralOrder);
            otherwise
                error('Unknown element type')
        end
        FEcell.Category{fi,1}=Econstant(ti).category;
        %FEcell.Type;
        FEcell.R_matrix{fi,1}=Element.R_matrix;
        FEcell.B_matrix{fi,1}=Element.B_matrix;
        FEcell.DetJac{fi,1}=Element.DetJac;
        FEcell.Nodes{fi,1}=Element.Nodes;
        FEcell.Elements{fi,1}=Element.Elements;
        FEcell.N_matrix{fi,1}=Element.N_matrix;
        %FEcell.Element_id;
        %FEcell.Material_id;
        FEcell.Field_dof_num{fi,1}=Econstant(ti).dof_num(field_dof_index);
    end
end
%==========================================================================
time=toc;fprintf("Finite element matrix is updated:%fs\n",time);
end