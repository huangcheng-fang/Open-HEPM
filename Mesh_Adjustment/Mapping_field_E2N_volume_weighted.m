%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function NEfield=Mapping_field_E2N_volume_weighted(Efield,Mesh)
Econstant=Element_constant_list();
EGauss=Element_Gauss_value();
Estress=Efield.Stress;
Estrain=Efield.Strain;
Epstrain=Efield.Pstrain;
Eid=Efield.Element_id;
nodes=Mesh.nodes;
elements=Mesh.elements;
%==========================================================================
container=Create_parallel_container(numel(Econstant));
for ti=1:1:numel(Econstant)
    eid=Eid{ti};
    if isempty(eid);continue;end
    node_num=Econstant(ti).node_num;
    var_num=Econstant(ti).var_num;
    pinvN=pinv(EGauss(ti).Interpolation);
    dN=EGauss(ti).Interpolation_derivative;
    weight=permute(EGauss(ti).Gauss_weight,[3,2,1,4]);
 
    enid=elements(eid,1:node_num).';
    element_nodes=pagetranspose(reshape(nodes(enid,:).',3,node_num,1,[]));
    Jac=pagemtimes(dN,element_nodes);
    detJac=PageDet(Jac);
    detJac=pagemtimes(detJac,weight);
    volume=sum(detJac,3);volume(:)=1;

    I=repmat(permute(enid,[1,4,3,2]),1,var_num+1);
    J=repmat(1:var_num+1,node_num,1,1,numel(eid));

    value=pagemtimes(pinvN,permute(Estress{ti},[3,1,2,4]));
    Nstress=[repmat(volume,node_num,1),pagemtimes(value,volume)];
 
    
    value=pagemtimes(pinvN,permute(Estrain{ti},[3,1,2,4]));
    Nstrain=[repmat(volume,node_num,1),pagemtimes(value,volume)];

    value=pagemtimes(pinvN,permute(Epstrain{ti},[3,1,2,4]));
    Npstrain=[repmat(volume,node_num,1),pagemtimes(value,volume)];

    container(ti,1).c1=[I(:),J(:)];
    container(ti,1).c2=Nstress(:);
    container(ti,1).c3=Nstrain(:);
    container(ti,1).c4=Npstrain(:);
end
[I,V1,V2,V3]=Merge_container(container,1);
Nstress=accumarray(I,V1);
Nstrain=accumarray(I,V2);
Npstrain=accumarray(I,V3);
%===========================Output=========================================
NEfield.Stress=(Nstress(:,1).\Nstress(:,2:end)).';
NEfield.Strain=(Nstrain(:,1).\Nstrain(:,2:end)).';
NEfield.Pstrain=(Npstrain(:,1).\Npstrain(:,2:end)).';
end