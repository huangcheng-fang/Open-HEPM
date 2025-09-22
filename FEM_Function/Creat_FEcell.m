%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FEcell=Creat_FEcell(Etype,MID,maxnum)
coder.varsize('Type')
Type=coder.nullcopy(cell(1,1));
coder.varsize('Element_id')
Element_id=coder.nullcopy(cell(1,1));
coder.varsize('Material_id')
Material_id=coder.nullcopy(cell(1,1));
Etype_list=unique(Etype);
MID_list=unique(MID);
%==========================================================================
ii=1;
for fi=1:1:numel(Etype_list)
    for mi=1:1:numel(MID_list)
        eid=find(Etype==Etype_list(fi)&MID==MID_list(mi));
        num=numel(eid);
        if isempty(eid);continue;end
        for k=1:1:ceil(num/maxnum)
            flag1=(k-1)*maxnum+1;
            flag2=min(flag1+maxnum-1,num);
            Element_id{ii,1}=eid(flag1:flag2);
            Type{ii,1}=Etype_list(fi);
            Material_id{ii,1}=MID_list(mi);
            ii=ii+1;
        end
    end
end
%==========================================================================
blank=coder.nullcopy(cell(numel(Element_id),1));
FEcell.Category=blank;
FEcell.Type=Type;
FEcell.R_matrix=blank;
FEcell.B_matrix=blank;
FEcell.DetJac=blank;
FEcell.Nodes=blank;
FEcell.Elements=blank;
FEcell.N_matrix=blank;
FEcell.Element_id=Element_id;
FEcell.Material_id=Material_id;
FEcell.Field_dof_num=blank;
end