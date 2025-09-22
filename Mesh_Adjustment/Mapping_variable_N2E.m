%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Result=Mapping_variable_N2E(NResult,Result,Mesh)
elements=Mesh.elements;
Etype=Mesh.etype;
Enum=size(Mesh.elements,1);
EIPnum=size(Result.Estress,2);
%=============================Inintial=====================================
Result.Estress=zeros(6,1,EIPnum,Enum);
Result.Estrain=zeros(6,1,EIPnum,Enum);
Result.Epstrain=zeros(6,1,EIPnum,Enum);
%=========================4-node tetrahedron element=======================
tet4loc=Etype==334;tet4e=elements(tet4loc,1:4);
Result.Estress(:,1,1,tet4loc)=Get_tetrahedron_elemental_result(NResult.Nstress,tet4e);
Result.Estrain(:,1,1,tet4loc)=Get_tetrahedron_elemental_result(NResult.Nstrain,tet4e);
Result.Epstrain(:,1,tet4loc)=Get_tetrahedron_elemental_result(NResult.Npstrain,tet4e);
end

function Eresult=Get_tetrahedron_elemental_result(Nresult,elements)
Eresult=0.25*(Nresult(:,elements(:,1))+Nresult(:,elements(:,2))+Nresult(:,elements(:,3))+Nresult(:,elements(:,4)));
end