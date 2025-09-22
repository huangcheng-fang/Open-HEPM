%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [elementCell,etypeCell,nodeCell,partIdCell]=Extract_mesh(FEMesh)
num=numel(FEMesh);
elementCell=coder.nullcopy(cell(num,1));
etypeCell=coder.nullcopy(cell(num,1));
nodeCell=coder.nullcopy(cell(num,1));
partIdCell=coder.nullcopy(cell(num,1));
for fi=1:1:num
    if strcmp(FEMesh(fi).GeometricOrder,'linear')
        type=334;
    else
        type=3310;
    end
    elementCell{fi}=FEMesh(fi).Elements.';
    etypeCell{fi}=repmat(type,size(elementCell{fi},1),1);
    nodeCell{fi}=FEMesh(fi).Nodes.';
    partIdCell{fi}=fi;
end
end