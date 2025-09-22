%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function B=Merge_cell(A,id)
index=zeros(numel(id)+1,1);
for ci=1:1:numel(id)
    index(ci+1)=index(ci)+numel(A{id(ci)});
end
B=zeros(index(end),1);
for ci=1:1:numel(id)
    B(index(ci)+1:index(ci+1))=A{id(ci)}(:);
end
B=unique(B);
end