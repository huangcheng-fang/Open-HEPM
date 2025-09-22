%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [C1,C2,C3,C4,C5,C6]=Merge_container(container,merge_dim)
if isempty(container)
    C1=[];C2=[];C3=[];C4=[];C5=[];C6=[];
    return
end
switch merge_dim
    case 1
        for ci=1:1:numel(container)
            container(ci).c1=container(ci).c1.';
            container(ci).c2=container(ci).c2.';
            container(ci).c3=container(ci).c3.';
            container(ci).c4=container(ci).c4.';
            container(ci).c5=container(ci).c5.';
            container(ci).c6=container(ci).c6.';
        end
    case 2
        
    otherwise
        error('Unsupported merge dimension')
end
index=zeros(numel(container)+1,6);
for ci=1:1:numel(container)
    index(ci+1,1)=index(ci,1)+size(container(ci).c1,2);
    index(ci+1,2)=index(ci,2)+size(container(ci).c2,2);
    index(ci+1,3)=index(ci,3)+size(container(ci).c3,2);
    index(ci+1,4)=index(ci,4)+size(container(ci).c4,2);
    index(ci+1,5)=index(ci,5)+size(container(ci).c5,2);
    index(ci+1,6)=index(ci,6)+size(container(ci).c6,2);
end
C1=zeros(size(container(1).c1,1),index(end,1));
C2=zeros(size(container(1).c2,1),index(end,2));
C3=zeros(size(container(1).c3,1),index(end,3));
C4=zeros(size(container(1).c4,1),index(end,4));
C5=zeros(size(container(1).c5,1),index(end,5));
C6=zeros(size(container(1).c6,1),index(end,6));
for ci=1:1:numel(container)
    C1(1:size(container(ci).c1,1),index(ci,1)+1:index(ci+1,1))=container(ci).c1;
    C2(1:size(container(ci).c2,1),index(ci,2)+1:index(ci+1,2))=container(ci).c2;
    C3(1:size(container(ci).c3,1),index(ci,3)+1:index(ci+1,3))=container(ci).c3;
    C4(1:size(container(ci).c4,1),index(ci,4)+1:index(ci+1,4))=container(ci).c4;
    C5(1:size(container(ci).c5,1),index(ci,5)+1:index(ci+1,5))=container(ci).c5;
    C6(1:size(container(ci).c6,1),index(ci,6)+1:index(ci+1,6))=container(ci).c6;
end
if merge_dim==1
C1=C1.';C2=C2.';C3=C3.';
C4=C4.';C5=C5.';C6=C6.';
end
end