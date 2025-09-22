%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function ips=Quadrilateral_inverse_interpolation(element_nodes,points)
ips=zeros(size(points,1),2);
for i=1:1:size(points,1)
    coor=points(i,:);
    ip = zeros(1,2);
    iter  = 100;

    inc = 1;
    while (inc < iter)
        N=Quadrilateral_interpolation(ip);
        dN=Quadrilateral_interpolation_derivative(ip);
        J=dN*element_nodes;
        f = N*element_nodes - coor;

        d_ip=f/J;
        ip=ip-d_ip;

        inc  = inc + 1;
        if norm(d_ip)<1e-10
            break
        end
    end
    if inc>=iter-2
        error('Exceeded the maximum number of iterations in calculating local coordinates')
    end
    ips(i,:)=ip;
end
end