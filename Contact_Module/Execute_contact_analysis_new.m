%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Contact_pair,ContactEquation]=Execute_contact_analysis_new(Contact_pair,ContactEquation,Contact_condition,Mesh,Set,stepU,locking)
maxdof=numel(stepU);positions=Mesh.positions;
%==========================================================================
Contact_pair=Judge_contact_pair_update(Contact_pair,Contact_condition,positions,1e-5);
Contact_pair=Update_contact_pair(Contact_pair,Contact_condition,Mesh,Set,stepU);
Contact_pair=Update_contact_field_new(Contact_pair,Contact_condition,positions,stepU,locking);
ContactEquation=Update_contact_Equation_new(ContactEquation,Contact_pair,Contact_condition,maxdof);
end