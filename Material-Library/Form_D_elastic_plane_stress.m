%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function D_elastic=Form_D_elastic_plane_stress(E,mu)
%简易版，还需要修改
A=E/(1-mu^2);B=A*mu;
D_elastic=[A,B,0,0;
           B,A,0,0;
           0,0,0,0;
           0,0,0,0.5*A*(1-mu)];
end