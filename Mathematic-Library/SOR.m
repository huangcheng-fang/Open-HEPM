%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [y,flag,err,iternum] = SOR(A,b,w,tol,maxnum)
% input: A 的对角线元素均不为 0 r: 松弛因子 e: 精度 M: 最大计算次数
% output: y:方程的解
n = length(A);
x0 = zeros(n,1);
y = zeros(n,1);

D = diag(diag(A),0);
U = triu(-A,1); %上三角矩阵
L = tril(-A,-1); %下三角矩阵
B = (D-w*L)\((1-w)*D+w*U);
f = w*((D-w*L)\b);
err=1;flag=0;iternum=0;
while err>tol
    y = B*x0+f;
    err=sum(abs(y-x0))/sum(abs(y));
    x0= y;iternum=iternum+1;
    if iternum>maxnum
        flag=1;
       break;
    end
end
end