%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FEcell=FEM_Element(Mesh,field_dof_index)
tic;
nodes=Mesh.nodes;
elements=Mesh.elements;
Etype=Mesh.etype;
Eactivation=Mesh.eactivation;
Econstant=Element_constant_list();
%==========================================================================
tnum=numel(Econstant);
blank2d=cell(tnum,1);
blank3d=cell(tnum,1);
blank4d=cell(tnum,1);
blankchar=cell(tnum,1);
coder.varsize('blankv2d')
coder.varsize('blankv3d')
coder.varsize('blankv4d')
coder.varsize('blankchar')
blankv2d=zeros(tnum-tnum,tnum-tnum);
blankv3d=zeros(tnum-tnum,tnum-tnum,tnum-tnum);
blankv4d=zeros(tnum-tnum,tnum-tnum,tnum-tnum,tnum-tnum);
blankvchar='';
for i=1:1:tnum;blank2d{i}=blankv2d;blank3d{i}=blankv3d;blank4d{i}=blankv4d;blankchar{i}=blankvchar;end
FEstruct=struct('Category',blankchar,'Type',blank2d,'R_matrix',blank4d,'B_matrix',blank4d,'DetJac',blank4d,'Nodes',blank4d,'Elements',blank4d,'N_matrix',blank3d,'Element_id',blank2d,'Material_id',blank2d,'Field_dof_num',blank2d);
%==============================main function===============================
active_type=false(tnum,1);
for ti=1:1:tnum
    loc=find(Etype==Econstant(ti).type);
    loc=loc(Eactivation(loc));
    if ~isempty(loc)
        active_type(ti)=true;
        Enid=elements(loc,1:Econstant(ti).node_num);
        IntegralOrder=Econstant(ti).int_order;
        switch Econstant(ti).type
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
        Element.Category=Econstant(ti).category;
        Element.Type=Econstant(ti).type;
        Element.Element_id=loc;
        Element.Field_dof_num=Econstant(ti).dof_num(field_dof_index);
        Element.Material_id=nan;
        FEstruct(ti,1)=Element;
    end
end
%==========================================================================
temp=struct2cell(FEstruct);
FEcell.Category=temp(1,active_type)';
FEcell.Type=temp(2,active_type)';
FEcell.R_matrix=temp(3,active_type)';
FEcell.B_matrix=temp(4,active_type)';
FEcell.DetJac=temp(5,active_type)';
FEcell.Nodes=temp(6,active_type)';
FEcell.Elements=temp(7,active_type)';
FEcell.N_matrix=temp(8,active_type)';
FEcell.Element_id=temp(9,active_type)';
FEcell.Material_id=temp(10,active_type)';
FEcell.Field_dof_num=temp(11,active_type)';
%==========================================================================
time=toc;fprintf("Finite element matrix is updated:%fs\n",time);
end