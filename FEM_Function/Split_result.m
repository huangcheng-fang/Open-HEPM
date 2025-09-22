%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Result=Split_result(Result,X,maxdof)
Result.dU=X(1:maxdof,1);
Result.U=Result.U+Result.dU;
Result.stepU=Result.stepU+Result.dU;
try
    dP=X(maxdof+1:maxdof+maxdof/3,1);
    Result.dP=dP;
    Result.P=Result.P+dP;
catch
end
end