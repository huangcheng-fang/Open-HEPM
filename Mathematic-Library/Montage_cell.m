%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function newData=Montage_cell(dataCell,dim)
if issparse(dataCell{end})
   error('Montage_cell does not surpport the sparse matrix')
end
if dim==2
    for j=1:1:numel(dataCell)
        dataCell{j}=dataCell{j}.';
    end
elseif dim~=1
        error('Montage_cell only surpports the 2D matrix')
end
%--------------------------------------------------------------------------
num=0;
maxsize2=0;
for j=1:1:numel(dataCell)
    maxsize2=max(maxsize2,size(dataCell{j},2));
    num=num+size(dataCell{j},1);
end
if ~islogical(dataCell{end})
    newData=zeros(num,maxsize2);
else
    newData=false(num,maxsize2);
end
%--------------------------------------------------------------------------
flag=0;
for j=1:1:numel(dataCell)
    newData(flag+1:flag+size(dataCell{j},1),1:size(dataCell{j},2))=dataCell{j};
    flag=flag+size(dataCell{j},1);
end
%--------------------------------------------------------------------------
if dim==2
    newData=newData.';
end
end