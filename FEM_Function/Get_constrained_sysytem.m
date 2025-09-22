%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [K,F]=Get_constrained_sysytem(K,F,constraint,CE)
tic
%--------------------------------------------------------------------------
K=K+CE.Penalty.matrix+CE.interface_storage;
F=F-CE.Penalty.vector-CE.interface_storage_vector;
%--------------------------------------------------------------------------
K=CE.Dual_Lagrange.matrix.'*K*CE.Dual_Lagrange.matrix+diag(CE.Dual_Lagrange.perturbation);
F=CE.Dual_Lagrange.matrix.'*(F-CE.Dual_Lagrange.Fcon);
loc=ismember(constraint(:,1),CE.Dual_Lagrange.slave_dof);
constraint=[constraint(~loc,:);CE.Dual_Lagrange.constraint];
%--------------------------------------------------------------------------
K=[K,CE.Lagrange.';CE.Lagrange,CE.perturbation];
F=[F;full(CE.perturbation_vector)];
%--------------------------------------------------------------------------
% contact_nid=unique(floor((Contact_matrix.constraint(:,1)+2)/3));
% loc=ismember(constraint(:,1),[contact_nid*3;contact_nid*3-1;contact_nid*3-2]);
% constraint=[constraint(~loc,:);Contact_matrix.constraint];
%------------------------------------------------------------------------------
PU=zeros(size(K,1),1);
PU(constraint(:,1),1)=constraint(:,2);
F=F-K*PU;F(constraint(:,1))=constraint(:,2);
Aeye=sparse(constraint(:,1),constraint(:,1),1,numel(F),numel(F));
Ceye=speye(numel(F))-Aeye;
K=Ceye*K*Ceye;K=K+Aeye;
%------------------------------------------------------------------------------
% loc=find(diag(K)==0);
% Aeye=sparse(loc,loc,1,numel(F),numel(F));
% F(loc)=0;K=K+Aeye;
time=toc;fprintf("Constrained_sysytem is calculated: %fs\n",time);
end