%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function local_coor=Triangular_inverse_interpolation(element_node,point)
k=[element_node(2,1)-element_node(1,1), element_node(3,1)-element_node(1,1);
    element_node(2,2)-element_node(1,2), element_node(3,2)-element_node(1,2);];
f=[point(:,1)'-element_node(1,1);point(:,2)'-element_node(1,2);];
local_coor=(k\f)';
end