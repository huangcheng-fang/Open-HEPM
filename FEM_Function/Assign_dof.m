%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FEdof=Assign_dof(CElements,CField_dof_num,DOF)
coder.varsize('FEdof')
FEdof=coder.nullcopy(cell(numel(CElements),numel(DOF)));
for ti=1:1:numel(CElements)
    Enid=CElements{ti};
    if isempty(Enid);continue;end
    for di=1:1:numel(DOF)
        dof_num=CField_dof_num{ti}(di);
        FEdof{ti,di}=reshape(DOF{di}(1:dof_num,Enid(:)),[],1,1,size(Enid,4));
    end
end
end