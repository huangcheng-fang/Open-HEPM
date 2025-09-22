%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function index=Parallel_data_index(Tnum,Cnum)
num=floor(Tnum/Cnum);
index=zeros(Cnum,2);
for ci=1:1:Cnum
    index(ci,:)=[(ci-1)*num+1,ci*num];
end
index(end)=max(index(end),Tnum);
end