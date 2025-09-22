%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [errorU,errorS]=Convergence_calculate_for_FEM(Mesh,Result,analyticU,analyticS)
field_dof_index=[1,2];
FEcell=Creat_FEcell(Mesh.etype,1,1e50);
FEcell=Calculate_FEcell(FEcell,Mesh,field_dof_index);
DOF=Get_node_dof_list(size(Mesh.nodes,1),{3});
FEdof=Assign_dof(FEcell.Elements,FEcell.Field_dof_num,DOF);
Efieldcell=Assign_field_result(Result.Efield,FEcell.Element_id);

N_matrix=FEcell.N_matrix;
Nodes=FEcell.Nodes;
DetJac=FEcell.DetJac;
U=Result.Nfield.U;

for fi=1:1:numel(N_matrix)
    ip_coor=pagemtimes(N_matrix{fi},Nodes{fi});
    x=ip_coor(:,1,:,:);y=ip_coor(:,2,:,:);z=ip_coor(:,3,:,:);
    ip_U_analytic=[];
    for ii=1:1:numel(analyticU)
        ip_U_analytic=[ip_U_analytic;eval(vectorize(analyticU{ii}))];
    end
    N_U=reshape(U(FEdof{fi}),size(Nodes{fi},2),size(Nodes{fi},1),1,size(Nodes{fi},4));
    N_U=pagetranspose(N_U);
    ip_U_numerical=pagetranspose(pagemtimes(N_matrix{fi},N_U));

    ip_S_analytic=[];
    for ii=1:1:numel(analyticS)
        ip_S_analytic=[ip_S_analytic;eval(vectorize(analyticS{ii}))];
    end
    ip_S_numerical=Efieldcell.Stress{fi};

    errorU=sum(pagemtimes((ip_U_analytic-ip_U_numerical).^2,DetJac{fi}),'all');
    errorU=errorU/sum(pagemtimes((ip_U_analytic).^2,DetJac{fi}),'all');
    errorU=sqrt(errorU);

    errorS=sum(pagemtimes((ip_S_analytic-ip_S_numerical).^2,DetJac{fi}),'all');
    errorS=errorS/sum(pagemtimes((ip_S_analytic).^2,DetJac{fi}),'all');
    errorS=sqrt(errorS);
end
end