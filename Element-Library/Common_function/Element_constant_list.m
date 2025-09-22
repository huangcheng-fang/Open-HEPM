%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%dof_num: The number of degrees of freedom at each node for different physical fields.
%         For instance,[3,1], 3 can represent the dof number of displacement, while 1 denotes that of pore pressure/temperature.
function [Econstant,Etype_list]=Element_constant_list()
Key=Keywords();solid=Key.solid;structure=Key.structure;
%==========================================================================
category  =     {solid; solid; solid; solid; structure; structure; structure};
type      =     {334; 338; 3381; 3310; 312; 323; 324};
node_num  =     {4; 8; 8; 10; 2; 3; 4};
dof_num   =     {[3,1]; [3,1]; [3,1]; [3,1]; [3,1]; [3,1]; [3,1]};
var_num   =     {6; 6; 6; 6; 1; 3; 3};
int_order =     {1; 2; 2; 2; 1; 1; 2};
int_point_num = {1; 8; 8; 4; 2; 1; 4};
%=============================Output=======================================
Econstant =struct('category',category,'type',type,'node_num',node_num,'dof_num',dof_num,'var_num',var_num,'int_order',int_order,'int_point_num',int_point_num);
Etype_list=zeros(numel(type),1);for i=1:1:numel(Etype_list);Etype_list(i)=type{i};end
end