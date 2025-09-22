%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function N=Triangular_quadratic_interpolation(ip)
x=ip(:,1);y=ip(:,2);
N=[2*(0.5-x-y).*(1-x-y), 2*x.*(x-0.5), 2*y.*(y-0.5), 4*x.*(1-x-y), 4*x.*y, 4*y.*(1-x-y)];
end