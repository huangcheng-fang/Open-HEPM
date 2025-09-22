%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function container=Create_parallel_container(num)
blank=cell(num,1);
for i=1:1:num
    blank{i}=zeros(num-num);
end
container=struct('c1',blank,'c2',blank,'c3',blank,'c4',blank,'c5',blank,'c6',blank);
end