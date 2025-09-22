%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function newDataCell=Collect_cell4D(dataCell,id,id_list)
csize=size(dataCell{1},[1,2,3]);
newDataCell=coder.nullcopy(cell(numel(id_list),1));
for i=1:1:numel(id_list)
    num=0;
    for j=1:1:numel(dataCell)
        if id{j}==id_list{i}
           num=num+size(dataCell{j},1);
        end
    end
    newDataCell{i}=zeros(csize(1),csize(2),csize(3),num);
end
%--------------------------------------------------------------------------
for i=1:1:numel(id_list)
    flag=0;
    for j=1:1:numel(dataCell)
        if id{j}==id_list{i}
            newDataCell{i}(:,:,:,flag+1:flag+size(dataCell{j},4))=dataCell{j};
            flag=flag+size(dataCell{j},4);
        end
    end
end
end