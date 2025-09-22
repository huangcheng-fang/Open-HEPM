%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Contact_pair,ContactEquation]=Initial_contact_solver(num,NCstress,NCstatus)
blank=cell(num,1);sparseblank=cell(num,1);logicblank=cell(num,1);logicvalue=cell(num,1);
for ci=1:1:num
    blank{ci}=zeros(num-num,1);
    logicblank{ci}=true(num-num,1);
    logicvalue{ci}=true;
    sparseblank{ci}=sparse(num-num,num-num);
end
Contact_pair=struct('area_vector',blank,'slave_dof',blank,'contact_matrix',sparseblank,'R_matrix',sparseblank,'contact_gap',blank,'contact_gap_P',blank,'contact_stress',blank,'contact_flux',blank,'contact_status',logicblank,'NCstress',blank,'NCstatus',logicblank,'update',logicvalue,'status_change',logicvalue);
%==========================================================================
for ci=1:1:num
    Contact_pair(ci).NCstress=NCstress(:,ci);
    Contact_pair(ci).NCstatus=NCstatus(:,ci);
end
%==========================================================================
maxdof=size(NCstress,1);
Penalty.matrix=sparse(maxdof,maxdof);
Penalty.vector=zeros(maxdof,1);
Dual_Lagrange.matrix=speye(maxdof);
Dual_Lagrange.Fcon=zeros(maxdof,1);
Dual_Lagrange.slave_dof=zeros(0,1);
Dual_Lagrange.constraint=zeros(0,2);
Dual_Lagrange.perturbation=sparse(maxdof,1);
%--------------------------------------------------------------------------
ContactEquation.Penalty=Penalty;
ContactEquation.Dual_Lagrange=Dual_Lagrange;
ContactEquation.update=true;
end