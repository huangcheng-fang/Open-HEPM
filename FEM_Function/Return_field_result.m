%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Efield=Return_field_result(Efield,Efieldcell)
id_index=Efield.id_index;
CEid=Efieldcell.Element_id;
for fi=1:1:numel(CEid)
    eid=CEid{fi};
    if isempty(eid);continue;end
    loc1=id_index(eid(1),1);
    loc2=id_index(eid,2);
    Efield.Stress{loc1}(:,:,:,loc2)=Efieldcell.Stress{fi};
    Efield.Strain{loc1}(:,:,:,loc2)=Efieldcell.Strain{fi};
    Efield.Pstrain{loc1}(:,:,:,loc2)=Efieldcell.Pstrain{fi};
    sz=size(Efieldcell.Parameter{fi},[1,2]);
    Efield.Parameter{loc1}(1:sz(1),1:sz(2),:,loc2)=Efieldcell.Parameter{fi};
end
end