%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function x=Truncated_errors(x,er)
n=floor(log10(abs(x))-er);
n(n==-inf)=0;
x=round(x./10.^n);
x=x.*10.^n;
x(abs(x)<1e-12)=0;
x=floor(x*1e12)/1e12;
end