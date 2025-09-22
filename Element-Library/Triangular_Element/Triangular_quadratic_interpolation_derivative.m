%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function dN=Triangular_quadratic_interpolation_derivative(ip)
ip=permute(ip,[2,3,1]);
x=4*ip(1,:,:); y=4*ip(2,:,:); 
OO=zeros(1,1,size(x,3));
%--------------------------------------------------------------------------
dN=[x + y - 3, x - 1, OO, 4 - y - 2*x, y, -y;
    x + y - 3, OO, y - 1, -x, x, 4 - 2*y - x];
end