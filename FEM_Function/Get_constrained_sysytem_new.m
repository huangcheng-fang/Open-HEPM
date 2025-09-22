%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [K,F]=Get_constrained_sysytem_new(K,F,constraint,U,CE,deactive_dof)
tic
constraint(:,2)=constraint(:,2)-U(constraint(:,1));
%--------------------------------------------------------------------------
K=K+CE.Penalty.matrix;
F=F-CE.Penalty.vector;
%--------------------------------------------------------------------------
diagK=diag(K);
penK=diagK(constraint(:,1))*1e8;
F(constraint(:,1))=F(constraint(:,1))+penK.*constraint(:,2);
Kpenalty=sparse(constraint(:,1),constraint(:,1),penK,numel(F),numel(F));
K=K+Kpenalty;
%--------------------------------------------------------------------------
K=CE.Dual_Lagrange.matrix.'*K*CE.Dual_Lagrange.matrix+CE.Dual_Lagrange.perturbation;
F=CE.Dual_Lagrange.matrix.'*(F-CE.Dual_Lagrange.Fcon);
% constraint=CE.Dual_Lagrange.constraint;
%--------------------------------------------------------------------------
% PU=zeros(size(K,1),1);
% PU(constraint(:,1),1)=constraint(:,2);
% F=F-K*PU;F(constraint(:,1))=constraint(:,2);
% Aeye=sparse(constraint(:,1),constraint(:,1),1,numel(F),numel(F));
% Ceye=speye(numel(F))-Aeye;
% K=Ceye*K*Ceye;K=K+Aeye;
%--------------------------------------------------------------------------
K_deactive=zeros(size(K,1),1);
K_deactive(deactive_dof)=1;
K_deactive=diag(sparse(K_deactive));
K=K+K_deactive;
%------------------------------------------------------------------------------
time=toc;fprintf("Constrained_sysytem is calculated: %fs\n",time);
end

% function [K,F]=Get_constrained_sysytem_new(K,F,constraint,U,CE,deactive_dof)
% tic
% constraint(:,2)=constraint(:,2)-U(constraint(:,1));
% %--------------------------------------------------------------------------
% K=K+CE.Penalty.matrix;
% F=F-CE.Penalty.vector;
% %--------------------------------------------------------------------------
% K=CE.Dual_Lagrange.matrix.'*K*CE.Dual_Lagrange.matrix+diag(CE.Dual_Lagrange.perturbation);
% F=CE.Dual_Lagrange.matrix.'*(F-CE.Dual_Lagrange.Fcon);
% loc=ismember(constraint(:,1),CE.Dual_Lagrange.slave_dof);
% constraint=[constraint(~loc,:);CE.Dual_Lagrange.constraint];
% %------------------------------------------------------------------------------
% PU=zeros(size(K,1),1);
% PU(constraint(:,1),1)=constraint(:,2);
% F=F-K*PU;F(constraint(:,1))=constraint(:,2);
% Aeye=sparse(constraint(:,1),constraint(:,1),1,numel(F),numel(F));
% Ceye=speye(numel(F))-Aeye;
% K=Ceye*K*Ceye;K=K+Aeye;
% %--------------------------------------------------------------------------
% K_deactive=zeros(size(K,1),1);
% K_deactive(deactive_dof)=1;
% K_deactive=diag(sparse(K_deactive));
% K=K+K_deactive;
% %------------------------------------------------------------------------------
% time=toc;fprintf("Constrained_sysytem is calculated: %fs\n",time);
% end