%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Efield=Mapping_result(Mesh0,Efield0,Mesh,Efield)
nodes0=Mesh0.nodes;
elements0=Mesh0.elements;
Einpart0=Mesh0.einpart;
Element_id0=Efield0.Element_id;
Stress0=Montage_cell_4D_to_2D(Efield0.Stress);
Strain0=Montage_cell_4D_to_2D(Efield0.Strain);
Pstrain0=Montage_cell_4D_to_2D(Efield0.Pstrain);
Parameter0=Montage_cell_4D_to_2D(Efield0.Parameter);
[Gauss_points0,Ginpart0]=Get_gauss_point(nodes0,elements0,Element_id0,Einpart0);
%--------------------------------------------------------------------------
nodes=Mesh.nodes;
elements=Mesh.elements;
Einpart=Mesh.einpart;
Element_id=Efield.Element_id;
Stress=Montage_cell_4D_to_2D(Efield.Stress);
Strain=Montage_cell_4D_to_2D(Efield.Strain);
Pstrain=Montage_cell_4D_to_2D(Efield.Pstrain);
Parameter=Montage_cell_4D_to_2D(Efield.Parameter);
[Gauss_points,Ginpart]=Get_gauss_point(nodes,elements,Element_id,Einpart);
%===========================Mapping Variable===============================
list=unique(Einpart0);
for ii=1:1:numel(list)
loc0=find(Ginpart0==list(ii));
loc=find(Ginpart==list(ii));
KDTree=KDTreeSearcher(Gauss_points0(loc0,:));
[PID1,Distance] = knnsearch(KDTree,Gauss_points(loc,:),'k',2);
Distance(Distance<1e-16)=1e-16;
weight=1./Distance;weight=sum(weight,2).\weight;
PID2=repmat((1:1:size(PID1,1))',1,size(PID1,2));
interp_matrix=sparse(PID1,PID2,weight,size(loc0,1),size(PID1,1));
Stress(:,loc)=Stress0(:,loc0)*interp_matrix;
Strain(:,loc)=Strain0(:,loc0)*interp_matrix;
Pstrain(:,loc)=Pstrain0(:,loc0)*interp_matrix;
Parameter(:,loc)=Parameter0(:,loc0)*interp_matrix;
end
%===============================Output=====================================
Econstant=Element_constant_list();
for ii=1:1:numel(Element_id)
    enum=numel(Element_id{ii});
    ip_num=Econstant(ii).int_point_num;
    var_num=size(Efield.Stress{ii},1);
    Efield.Stress{ii}(1:var_num,1,1:ip_num,1:enum)=Stress(1:var_num,1:ip_num*enum);
    Efield.Strain{ii}(1:var_num,1,1:ip_num,1:enum)=Strain(1:var_num,1:ip_num*enum);
    Efield.Pstrain{ii}(1:var_num,1,1:ip_num,1:enum)=Pstrain(1:var_num,1:ip_num*enum);
    Efield.Parameter{ii}(:,1,1:ip_num,1:enum)=Parameter(:,1:ip_num*enum);
end
end



function [Gauss_points,Ginpart]=Get_gauss_point(nodes,elements,Element_id,Einpart)
EGauss=Element_Gauss_value();
[Econstant,Etype_list]=Element_constant_list();
coder.varsize('Gauss_points')
coder.varsize('Ginpart')
Gauss_points=coder.nullcopy(cell(0,numel(Etype_list)));
Ginpart=coder.nullcopy(cell(0,numel(Etype_list)));
for ii=1:1:numel(Etype_list)
    eid=Element_id{ii};
    if isempty(eid);continue;end
    enid=elements(eid,1:Econstant(ii).node_num);
    CN_matrix=EGauss(ii).Interpolation;
    CNodes=reshape(nodes(enid.',1:3).',3,size(enid,2),1,size(enid,1));
    Gauss_points{ii}=pagemtimes(CN_matrix,'none',CNodes,'transpose');
    Ginpart{ii}=repmat(permute(Einpart(eid),[4,2,3,1]),1,1,size(Gauss_points{ii},3),1);
end
Gauss_points=Montage_cell_4D_to_2D(Gauss_points).';
Ginpart=Montage_cell_4D_to_2D(Ginpart).';
end
