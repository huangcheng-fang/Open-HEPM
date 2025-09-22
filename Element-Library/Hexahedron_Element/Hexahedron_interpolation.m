%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function N=Hexahedron_interpolation(ip)
X=ip(:,1);Y=ip(:,2);Z=ip(:,3);
auxX=[-1,-1,-1,-1,1,1,1,1];
auxY=[1,-1,-1,1,1,-1,-1,1];
auxZ=[1,1,-1,-1,1,1,-1,-1];
N=(1+X*auxX).*(1+Y*auxY).*(1+Z*auxZ)*0.125;
end