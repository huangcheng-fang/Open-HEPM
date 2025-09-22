%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function data=Montage_cell_4D_to_2D(dataCell)
for ii=1:1:numel(dataCell)
sz=size(dataCell{ii},[1,2,3,4]);
dataCell{ii}=reshape(dataCell{ii},sz(1)*sz(2),sz(3)*sz(4));
end
data=Montage_cell(dataCell,2);
end