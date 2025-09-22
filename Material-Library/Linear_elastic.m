%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%elasPara:[Elastic modulus;Poisson's ratio]
%dstrain:6x1xIPNUM*ENUM vector,[Ex;Ey;Ez;Exy;Eyz;Ezx]
%stress:6x1xIPNUM*ENUM vector,[Sx;Sy;Sz;Sxy;Syz;Szx]
function [D,stress]=Linear_elastic(elasPara,dstrain,stress)
E=elasPara(1);mu=elasPara(2);
%==========================================================================
K=E/(3*(1-2*mu));
G=E/(2*(1+mu));
G2=2*G;H=K-G2/3;F=G2+H;
D=[F, H, H, 0, 0, 0;
           H, F, H, 0, 0, 0;
           H, H, F, 0, 0, 0;
           0, 0, 0, G, 0, 0;
           0, 0, 0, 0, G, 0;
           0, 0, 0, 0, 0, G;];
%==========================================================================
stress=stress+pagemtimes(D,dstrain);
end