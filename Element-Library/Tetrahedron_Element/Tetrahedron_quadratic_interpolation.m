%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function N=Tetrahedron_quadratic_interpolation(ip)
x=ip(:,1); 
y=ip(:,2); 
z=ip(:,3); 
T=(1-x - y - z);
N = [2.*-T.*(0.5 - T), x.*(2*x - 1), y.*(2*y - 1), z.*(2*z - 1), 4*x.*T, 4*x.*y, 4*y.*T, 4*z.*T, 4*x.*z, 4*y.*z];
end