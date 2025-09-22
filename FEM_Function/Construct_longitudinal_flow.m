%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [KHM,FHM]=Construct_longitudinal_flow(interface_storage,Interface_element,Contact_pair,Mesh,Sset,Result)
Key=Keywords();
nodes=Mesh.nodes;
surfaces=Mesh.surfaces;
maxdof=numel(nodes);
%=========================Interface_element_thickness======================
global thickness0
thickness=zeros(maxdof/3,1);
for ci=1:numel(Contact_pair)
    slave_nid=Contact_pair(ci).slave_dof(3:3:end)/3;
    thickness(slave_nid)=thickness(slave_nid)+Contact_pair(ci).contact_gap;
end
thickness=max(-thickness,1e-8);
if isempty(thickness0);thickness0=thickness;end
thickness=(0.95*thickness+0.05*thickness0);
thickness(thickness<1e-3)=thickness0(thickness<1e-3);
%==========================Interface_element_matrix========================
KHM=interface_storage;
for ii=1:1:numel(Interface_element)
    eid=Merge_cell(Sset,Interface_element(ii).set);
    elements=surfaces(eid,1:3);
    parameter=Interface_element(ii).parameter;
    ini_thickness=parameter.ini_thickness;
    element_thickness=reshape(mean(thickness(elements),2),1,1,1,[]);
    element_thickness=max(element_thickness,ini_thickness);
    %----------------------------------------------------------------------
    De=Form_D_elastic_plane_strain(parameter.elastic_modulus,parameter.poisson_ratio);
    EDe=repmat(De([1,2,4],[1,2,4]),1,1,1,size(elements,1));
    D_permeability=diag(repmat(parameter.permeability,2,1));
    switch parameter.flow_type
        case Key.interface_cubic
            ED_permeability=pagemtimes(D_permeability,element_thickness.^2);
        case Key.interface_darcy
            ED_permeability=repmat(D_permeability,1,1,1,size(elements,1));
        otherwise
            error('Unkown interface flow type')
    end
    %----------------------------------------------------------------------
    Tri3D=Triangular_Element3D(elements,nodes);
    Tri3D.detJac=pagemtimes(Tri3D.detJac,reshape(element_thickness,1,1,[]));
    KS=From_longitudinal_flow_solid_matrix(Tri3D,EDe,maxdof);
    [KF,KCP]=From_longitudinal_flow_fluid_matrix(Tri3D,ED_permeability,maxdof);
    KHM=KHM+[KS,-KCP*0;-KCP.'*0,-KF];
end
FHM=[KHM(1:maxdof,:)*[Result.U;Result.P];KHM(maxdof+1:end,:)*[Result.stepU;Result.P]];
end