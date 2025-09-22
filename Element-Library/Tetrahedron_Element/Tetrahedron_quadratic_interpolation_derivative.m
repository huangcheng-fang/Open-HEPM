%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function dN=Tetrahedron_quadratic_interpolation_derivative(ip)
ip=permute(ip,[2,3,1]);
x=4*ip(1,:,:); y=4*ip(2,:,:); z=4*ip(3,:,:);
OO=zeros(1,1,size(x,3));
%--------------------------------------------------------------------------
dN=[x + y + z - 3, x - 1, OO, OO, 4 - y - z - 2*x, y, -y, -z, z, OO;
    x + y + z - 3, OO, y - 1, OO, -x, x, 4 - 2*y - z - x, -z, OO, z;
    x + y + z - 3, OO, OO, z - 1, -x, OO, -y, 4 - y - 2*z - x, x, y];
end