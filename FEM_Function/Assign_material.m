%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FEcell=Assign_material(FE,MID)
coder.varsize('blank')
blank=coder.nullcopy(cell(1,1));
FEcell.Category=blank;
FEcell.Type=blank;
FEcell.R_matrix=blank;
FEcell.B_matrix=blank;
FEcell.DetJac=blank;
FEcell.Nodes=blank;
FEcell.Elements=blank;
FEcell.N_matrix=blank;
FEcell.Element_id=blank;
FEcell.Material_id=blank;
FEcell.Field_dof_num=blank;

ii=1;
for fi=1:1:numel(FE.Type)
    mis=MID(FE.Element_id{fi},1);
    mino=unique(mis);
    if isempty(mis);continue;end
    for mi=1:1:numel(mino)
        loc=mis==mino(mi);
        FEcell.Category{ii,1}=FE.Category{fi};
        FEcell.Type{ii,1}=FE.Type{fi};
        FEcell.R_matrix{ii,1}=FE.R_matrix{fi}(:,:,:,loc);
        FEcell.B_matrix{ii,1}=FE.B_matrix{fi}(:,:,:,loc);
        FEcell.DetJac{ii,1}=FE.DetJac{fi}(:,:,:,loc);
        FEcell.Nodes{ii,1}=FE.Nodes{fi}(:,:,:,loc);
        FEcell.Elements{ii,1}=FE.Elements{fi}(:,:,:,loc);
        FEcell.N_matrix{ii,1}=FE.N_matrix{fi};
        FEcell.Material_id{ii,1}=mino(mi);
        FEcell.Element_id{ii,1}=FE.Element_id{fi}(loc,1);
        FEcell.Field_dof_num{ii,1}=FE.Field_dof_num{fi};
        ii=ii+1;
    end
end
end