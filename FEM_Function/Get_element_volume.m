%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Evolume=Get_element_volume(FE,Enum)
Evolume=zeros(Enum,1);
for fi=1:1:numel(FE.Element_id)
    Evolume(FE.Element_id{fi})=reshape(sum(FE.DetJac{fi},3),[],1);
end
end