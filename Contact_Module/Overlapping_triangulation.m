%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [List,Vertices]=Overlapping_triangulation(Vertices)
if size(Vertices,1)==3
    List=[1,2,3];
    return
end
a=size(Vertices,1);
vmid=sum(Vertices,1)/a;
List=[(1:a)',[(2:a)';1],repmat(a+1,a,1)];
Vertices=[Vertices;vmid];
end