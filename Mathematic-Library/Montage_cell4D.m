%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function newData=Montage_cell4D(dataCell)
num=0;
maxsize1=0;
maxsize2=0;
maxsize3=0;
for j=1:1:numel(dataCell)
    maxsize1=max(maxsize1,size(dataCell{j},1));
    maxsize2=max(maxsize2,size(dataCell{j},2));
    maxsize3=max(maxsize3,size(dataCell{j},3));
    num=num+size(dataCell{j},4);
end
newData=zeros(maxsize1,maxsize2,maxsize3,num);
%--------------------------------------------------------------------------
flag=0;
for j=1:1:numel(dataCell)
    size1=size(dataCell{j},1);
    size2=size(dataCell{j},2);
    size3=size(dataCell{j},3);
    size4=size(dataCell{j},4);
    newData(1:size1,1:size2,1:size3,flag+1:flag+size4)=dataCell{j};
    flag=flag+size4;
end
end