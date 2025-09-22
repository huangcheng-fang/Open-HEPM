%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function NEfield=Mapping_field_E2N(Efield,Mesh)
Econstant=Element_constant_list();
EGauss=Element_Gauss_value();
Estress=Efield.Stress;
Estrain=Efield.Strain;
Epstrain=Efield.Pstrain;
Eparameter=Efield.Parameter;
Eid=Efield.Element_id;
elements=Mesh.elements;
%==========================================================================
container=Create_parallel_container(numel(Econstant));
for ti=1:1:numel(Econstant)
    eid=Eid{ti};
    if isempty(eid);continue;end
    node_num=Econstant(ti).node_num;
    var_num=Econstant(ti).var_num;
    pinvN=pinv(EGauss(ti).Interpolation);
    if Econstant(ti).type==3310
        node_num=4;
        ip=Hammer_integral(3,2);
        pinvN=inv(Tetrahedron_interpolation(ip));
    end
 
    enid=elements(eid,1:node_num).';

    I=repmat(permute(enid,[1,4,3,2]),1,var_num+1);
    J=repmat(1:var_num+1,node_num,1,1,numel(eid));

    value=pagemtimes(pinvN,permute(Estress{ti},[3,1,2,4]));
    Nstress=[ones(node_num,1,1,numel(eid)),value];
    
    value=pagemtimes(pinvN,permute(Estrain{ti},[3,1,2,4]));
    Nstrain=[ones(node_num,1,1,numel(eid)),value];

    value=pagemtimes(pinvN,permute(Epstrain{ti},[3,1,2,4]));
    Npstrain=[ones(node_num,1,1,numel(eid)),value];

    value=pagemtimes(pinvN,permute(Eparameter{ti},[3,1,2,4]));
    Nparameter=[ones(node_num,1,1,numel(eid)),value];
    pvar_num=size(Eparameter{ti},1);
    Ip=repmat(permute(enid,[1,4,3,2]),1,pvar_num+1);
    Jp=repmat(1:pvar_num+1,node_num,1,1,numel(eid));

    container(ti,1).c1=[I(:),J(:)];
    container(ti,1).c2=Nstress(:);
    container(ti,1).c3=Nstrain(:);
    container(ti,1).c4=Npstrain(:);
    container(ti,1).c5=[Ip(:),Jp(:)];
    container(ti,1).c6=Nparameter(:);
end
[I,V1,V2,V3,I4,V4]=Merge_container(container,1);
matsize=[size(Mesh.nodes,1),size(Nstrain,2)];
Nstress=accumarray(I,V1,matsize);
Nstrain=accumarray(I,V2,matsize);
Npstrain=accumarray(I,V3,matsize);
Nparameter=accumarray(I4,V4,[size(Mesh.nodes,1),size(Nparameter,2)]);
%===========================Output=========================================
NEfield.Stress=(Nstress(:,1).\Nstress(:,2:end)).';
NEfield.Strain=(Nstrain(:,1).\Nstrain(:,2:end)).';
NEfield.Pstrain=(Npstrain(:,1).\Npstrain(:,2:end)).';
NEfield.Parameter=(Nparameter(:,1).\Nparameter(:,2:end)).';
end