%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function ip=Inverse_interpolation(type,element_node,point)
switch type
    case {3,323,223}
        ip=Triangular_inverse_interpolation(element_node,point);
    case {4,324,224}
        ip=Quadrilateral_inverse_interpolation(element_node,point);
    case 334
        ip=Tetrahedron_inverse_interpolation(element_node,point);
    case 338
        ip=Hexahedron_inverse_interpolation(element_node,point);
    case {312}
        ip=Line_inverse_interpolation(element_node,point);
    otherwise
        error('Unkonwn element type')
end
end