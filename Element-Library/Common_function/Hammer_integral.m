%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [ip,wt]=Hammer_integral(dim,order)
if dim==3
    switch order
        case 1
            ip=[1/4,1/4,1/4];wt=1/6;
        case 2
            ip=[0.58541020,0.13819660,0.13819660;0.13819660,0.58541020,0.13819660;0.13819660,0.13819660,0.58541020;0.13819660,0.13819660,0.13819660];wt=zeros(4,1)+1/24;
        case 3
            ip=[1/4,1/4,1/4;1/2,1/6,1/6;1/6,1/2,1/6;1/6,1/6,1/2;1/6,1/6,1/6];wt=[-0.8;0.45;0.45;0.45;0.45]/6;
        otherwise
            ip=nan;wt=nan;
            error('Unkonwn order')
    end
elseif dim==2
    switch order
        case 1
            ip=[1/3,1/3];wt=0.5;
        case 2
            ip=[2/3,1/6;1/6,2/3;1/6,1/6];wt=zeros(3,1)+1/6;
        case 3
            ip=[1/3,1/3;0.6,0.2;0.2,0.6;0.2,0.2];wt=[-27/96;25/96;25/96;25/96];
        otherwise
            ip=nan;wt=nan;
            error('Unkonwn order')
    end
else
    error('Unkonwn Hammer integral dimension')
end
end