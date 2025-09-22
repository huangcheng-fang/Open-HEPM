%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function pwn=Start_Parallel_Processor(pwn)
if pwn<2
    pwn=1;
    return;
end
paral= gcp('nocreate');
if ~isempty(paral)
    workers=paral.NumWorkers;
    if pwn~=workers
        delete(paral);
        paral=parpool(pwn);
        pwn=paral.NumWorkers;
    end
else
    parpool(pwn);
end
end    