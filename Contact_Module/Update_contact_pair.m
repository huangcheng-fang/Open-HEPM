%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Contact_pair=Update_contact_pair(Contact_pair,Contact_condition,Mesh,Set,stepU)
starttime=tic;
Stype=Mesh.stype;
surfaces=Mesh.surfaces;
Etype=Mesh.etype;
elements=Mesh.elements;
nodes=Mesh.nodes;
positions=Mesh.positions;
Sset=Set.surface_set;
Eset=Set.element_set;
chageflag=false;
%==========================================================================
for ci=1:numel(Contact_pair)
    if ~Contact_pair(ci).update
        continue
    end
    if Contact_condition(ci).large_sliding
        coors=positions;
    else
        coors=nodes;
    end
    switch Contact_condition(ci).discretization
        case 'surface_to_surface'
            plane_orientation=Contact_condition(ci).plane_orientation;
            slave=Merge_cell(Sset,Contact_condition(ci).slave_set);
            master=Merge_cell(Sset,Contact_condition(ci).master_set);
            slave_surface=surfaces(slave,:);
            master_surface=surfaces(master,:);
            slave_type=Stype(slave,:);
            master_type=Stype(master,:);
            isdual=strcmp(Contact_condition(ci).enforcement,'dual_lagrange');
            [contact_matrix,R_matrix,area_vector,slave_dof]=Surface_to_surface(slave_surface,master_surface,slave_type,master_type,coors,plane_orientation,true);
        case 'surface_to_volume'
            slave=Merge_cell(Sset,Contact_condition(ci).slave_set);
            master=Merge_cell(Eset,Contact_condition(ci).master_set);
            slave_surface=surfaces(slave,:);
            master_element=elements(master,:);
            slave_type=Stype(slave,:);
            master_type=Etype(master,:);
            isdual=strcmp(Contact_condition(ci).enforcement,'dual_lagrange');
            [contact_matrix,R_matrix,area_vector,slave_dof]=Surface_to_volume(slave_surface,master_element,slave_type,master_type,coors,isdual);
        case 'line_to_volume'
            slave=Merge_cell(Eset,Contact_condition(ci).slave_set);
            master=Merge_cell(Eset,Contact_condition(ci).master_set);
            slave_element=elements(slave,:);
            master_element=elements(master,:);
            slave_type=Etype(slave,:);
            master_type=Etype(master,:);
            isdual=strcmp(Contact_condition(ci).enforcement,'dual_lagrange');
            [contact_matrix,R_matrix,area_vector,slave_dof]=Line_to_volume(slave_element,master_element,slave_type,master_type,coors,isdual);
        case 'node_to_surface'
            plane_orientation=Contact_condition(ci).plane_orientation;
            slave=Merge_cell(Sset,Contact_condition(ci).slave_set);
            master=Merge_cell(Sset,Contact_condition(ci).master_set);
            slave_surface=surfaces(slave,:);
            master_surface=surfaces(master,:);
            [contact_matrix,R_matrix,area_vector,slave_dof]=node_to_surface(slave_surface,master_surface,coors,plane_orientation);
        
        otherwise
            contact_matrix=sparse(num-num,num-num);R_matrix=contact_matrix;area_vector=zeros(0,1);slave_dof=zeros(0,1);
            error(['Unknow contact type:',Contact_condition(ci).type]);
    end
%     contact_gap=R_matrix*(contact_matrix*stepU);
%     contact_gap(3:3:end)=R_matrix(3:3:end,:)*(contact_matrix*reshape(positions',[],1));
    Contact_pair(ci).area_vector=area_vector(:);
    Contact_pair(ci).slave_dof=slave_dof;
%     Contact_pair(ci).contact_gap=contact_gap;
%     Contact_pair(ci).contact_stress=R_matrix*Contact_pair(ci).NCstress(slave_dof)./area_vector;
%     Contact_pair(ci).contact_status=Contact_pair(ci).NCstatus(slave_dof);
    Contact_pair(ci).contact_matrix=contact_matrix;
    Contact_pair(ci).R_matrix=R_matrix;
    chageflag=true;
end
%==========================================================================
if chageflag
time=toc(starttime);fprintf("Contact_pair is updated:%fs\n",time);
end
end