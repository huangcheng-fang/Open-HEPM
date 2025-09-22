%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function newDataCell=Collect_cell(dataCell,id,id_list,dim)
if issparse(dataCell{end})
   error('Collect_cell does not surpport the sparse matrix')
end
if dim==2
    for j=1:1:numel(dataCell)
        dataCell{j}=dataCell{j}.';
    end
elseif dim~=1
        error('Collect_cell only surpports the 2D matrix')
end
%--------------------------------------------------------------------------
size2=size(dataCell{1},2);
newDataCell=coder.nullcopy(cell(numel(id_list),1));
for i=1:1:numel(id_list)
    num=0;
    for j=1:1:numel(dataCell)
        if id{j}==id_list{i}
           num=num+size(dataCell{j},1);
        end
    end
    newDataCell{i}=zeros(num,size2);
end
%--------------------------------------------------------------------------
for i=1:1:numel(id_list)
    flag=0;
    for j=1:1:numel(dataCell)
        if id{j}==id_list{i}
            newDataCell{i}(flag+1:flag+size(dataCell{j},1),:)=dataCell{j};
            flag=flag+size(dataCell{j},1);
        end
    end
end
%--------------------------------------------------------------------------
if dim==2
    for j=1:1:numel(newDataCell)
        newDataCell{j}=newDataCell{j}.';
    end
end
end