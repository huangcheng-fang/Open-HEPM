%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Convert_element_type(Model,type1,type2)
if isnan(type2)
    loc=Model.Mesh.etype==type1;
    Model.Mesh.etype(loc,:)=[];
    Model.Mesh.elements(loc,:)=[];
    Model.Mesh.einpart(loc,:)=[];
    Model.Mesh.ecenter(loc,:)=[];
    Model.Mesh.evolume(loc,:)=[];
else
   error('Unsupported type')
end
end