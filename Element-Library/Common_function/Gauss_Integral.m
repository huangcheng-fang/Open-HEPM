%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [integral_points,weights]=Gauss_Integral(dim,order)
[ip,wt]=grule(order);
integral_points=zeros(order^dim,dim);weights=zeros(order^dim,1);
switch dim
    case 1
        integral_points=ip;
        weights=wt;
    case 2
        flag=1;
        for i=1:1:order
            for j=1:1:order
                integral_points(flag,1:2)=[ip(i),ip(j)];
                weights(flag,1)=wt(i)*wt(j);
                flag=flag+1;
            end
        end
    case 3
        flag=1;
        for i=1:1:order
            for j=1:1:order
                for k=1:1:order
                    integral_points(flag,1:3)=[ip(i),ip(j),ip(k)];
                    weights(flag,1)=wt(i)*wt(j)*wt(k);
                    flag=flag+1;
                end
            end
        end
    otherwise
        error('Unknown Gauss integral dimension')
end
end