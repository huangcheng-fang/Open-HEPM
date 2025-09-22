%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function N=Get_shape_function(type,ip_local)
switch type
    case {3,323,223}
        N=Triangular_interpolation(ip_local);
    case {4,324,224}
        N=Quadrilateral_interpolation(ip_local);
    case {334}
        N=Tetrahedron_interpolation(ip_local);
    case {338}
        N=Hexahedron_interpolation(ip_local);
    case {312}
        N=Line_interpolation(ip_local);
    otherwise
        error('Unknown element type')
end
end