%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function CP=Get_contact_stress_from_solid(CP,Result,Mesh)
NResult=Mapping_variable_E2N(Result,Mesh);
Nstress=NResult.Nstress;
Nstress(isnan(Nstress))=0;
for ci=1:1:numel(CP)
    R_matrix=CP(ci).R_matrix;
    contact_matrix=CP(ci).contact_matrix(3:3:end,3:3:end);

    Nstress_tensor=reshape((contact_matrix*reshape(Nstress,6,[]).').',6,1,[]);
    Nstress_tensor=reshape(Nstress_tensor([1,4,6,4,2,5,6,5,3],1,:),3,3,[]);
    T=repmat(speye(3,3),1,size(R_matrix,1)/3);
    R_tensor=reshape(full(T*R_matrix),3,3,[]);

    Cstress=pagemtimes(Nstress_tensor,pagetranspose(R_tensor(3,:,:)));
    contact_stress=pagemtimes(R_tensor,Cstress);

    CP(ci).contact_stress=-contact_stress(:);
end
end