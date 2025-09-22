%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function CEfield=Assign_field_result(Efield,CEid)
coder.varsize('blank')
blank=coder.nullcopy(cell(numel(CEid),1));
CEfield.Stress=blank;
CEfield.Strain=blank;
CEfield.Pstrain=blank;
CEfield.Parameter=blank;
CEfield.Element_id=blank;
%==========================================================================
id_index=Efield.id_index;
for fi=1:1:numel(CEid)
    eid=CEid{fi};
    if isempty(eid);continue;end
    loc1=id_index(eid(1),1);
    loc2=id_index(eid,2);
    CEfield.Stress{fi}=Efield.Stress{loc1}(:,:,:,loc2);
    CEfield.Strain{fi}=Efield.Strain{loc1}(:,:,:,loc2);
    CEfield.Pstrain{fi}=Efield.Pstrain{loc1}(:,:,:,loc2);
    CEfield.Parameter{fi}=Efield.Parameter{loc1}(:,:,:,loc2);
    CEfield.Element_id{fi}=eid;
end
end