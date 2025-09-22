%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Pfield=Return_field_result_FPEM(CPfield)
Pfield.Stress=Montage_cell4D(CPfield.Stress);
Pfield.Strain=Montage_cell4D(CPfield.Strain);
Pfield.Pstrain=Montage_cell4D(CPfield.Pstrain);
Pfield.Parameter=Montage_cell4D(CPfield.Parameter);
end