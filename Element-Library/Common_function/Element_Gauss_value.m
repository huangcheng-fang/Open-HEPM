%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function EGauss=Element_Gauss_value()
Econstant = Element_constant_list();

[ip1, wt1] = Hammer_integral (3,Econstant(1).int_order);
[ip2, wt2] = Gauss_Integral  (3,Econstant(2).int_order);
[ip3, wt3] = Gauss_Integral  (3,Econstant(3).int_order);
[ip4, wt4] = Hammer_integral (3,Econstant(4).int_order);
[ip5, wt5] = Gauss_Integral  (1,Econstant(5).int_order);
[ip6, wt6] = Hammer_integral (2,Econstant(6).int_order);
[ip7, wt7] = Gauss_Integral  (2,Econstant(7).int_order);

Gauss_point   =  {ip1; ip2; ip3; ip4; ip5; ip6; ip7};

Gauss_weight  =  {wt1; wt2; wt3; wt4; wt5; wt6; wt7};

Interpolation =  {Tetrahedron_interpolation(ip1);
                  Hexahedron_interpolation(ip2);
                  Hexahedron_interpolation(ip3);
                  Tetrahedron_quadratic_interpolation(ip4);
                  Line_interpolation(ip5);
                  Triangular_interpolation(ip6);
                  Quadrilateral_interpolation(ip7);};

Interpolation_derivative = {Tetrahedron_interpolation_derivative(ip1);
                            Hexahedron_interpolation_derivative(ip2);
                            Hexahedron_interpolation_derivative(ip3);
                            Tetrahedron_quadratic_interpolation_derivative(ip4);
                            Line_interpolation_derivative(ip5);
                            Triangular_interpolation_derivative(ip6);
                            Quadrilateral_interpolation_derivative(ip7);};
%=============================Output=======================================
EGauss =struct('Gauss_point',Gauss_point,'Gauss_weight',Gauss_weight,'Interpolation',Interpolation,'Interpolation_derivative',Interpolation_derivative);
end