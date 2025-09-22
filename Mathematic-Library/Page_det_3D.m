%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function A=Page_det_3D(J)
A=J(1,1,:).*J(2,2,:).*J(3,3,:) - J(1,1,:).*J(3,2,:).*J(2,3,:) - J(2,1,:).*J(1,2,:).*J(3,3,:) + J(2,1,:).*J(3,2,:).*J(1,3,:) + J(3,1,:).*J(1,2,:).*J(2,3,:) - J(3,1,:).*J(2,2,:).*J(1,3,:);
A=A(:);
end