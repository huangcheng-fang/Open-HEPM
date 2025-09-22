%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Contact_pair,ContactEquation]=Execute_contact_analysis(Contact_pair,Contact_condition,Mesh,Result,locking,HMC)
NCstress=Result.NCstress;NCflux=Result.NCflux;dU=Result.dU;P=Result.P;U=Result.U;stepU=Result.stepU;maxdof=numel(dU);
%---------------------------------------------------------------------------------------------
Contact_pair=Update_contact_field(Contact_pair,Contact_condition,Mesh.positions,NCstress,NCflux,U,stepU,P,locking);
ContactEquation=Update_contact_Equation(Contact_pair,Contact_condition,maxdof);
if ~HMC;ContactEquation.interface_storage_vector=0;ContactEquation.interface_storage=sparse(0);return;end
FluidContactEquation=Update_fluid_contact_Equation(Contact_pair,Contact_condition,P,dU,maxdof/3);
%--------------------------------------------------------------------------
ContactEquation.Penalty.vector=[ContactEquation.Penalty.vector;FluidContactEquation.Penalty.vector];
ContactEquation.Penalty.matrix=[ContactEquation.Penalty.matrix,sparse(maxdof,maxdof/3);sparse(maxdof/3,maxdof),FluidContactEquation.Penalty.matrix];
%--------------------------------------------------------------------------
ContactEquation.Dual_Lagrange.matrix=[ContactEquation.Dual_Lagrange.matrix,sparse(maxdof,maxdof/3);sparse(maxdof/3,maxdof),FluidContactEquation.Dual_Lagrange.matrix];
ContactEquation.Dual_Lagrange.Fcon=[ContactEquation.Dual_Lagrange.Fcon;FluidContactEquation.Dual_Lagrange.Fcon];
ContactEquation.Dual_Lagrange.slave_dof=[ContactEquation.Dual_Lagrange.slave_dof;FluidContactEquation.Dual_Lagrange.slave_dof];
ContactEquation.Dual_Lagrange.constraint=[ContactEquation.Dual_Lagrange.constraint;FluidContactEquation.Dual_Lagrange.constraint];
ContactEquation.Dual_Lagrange.perturbation=[ContactEquation.Dual_Lagrange.perturbation;FluidContactEquation.Dual_Lagrange.perturbation];
%--------------------------------------------------------------------------
ContactEquation.Fcon=[ContactEquation.Fcon;FluidContactEquation.Fcon];
ContactEquation.constraint=[ContactEquation.constraint;FluidContactEquation.constraint];
temp1=sparse(size(ContactEquation.Lagrange,1),maxdof/3);
temp2=sparse(size(FluidContactEquation.Lagrange,1),maxdof);
ContactEquation.Lagrange=[ContactEquation.Lagrange,temp1;temp2,FluidContactEquation.Lagrange];
ContactEquation.perturbation=diag(sparse([ContactEquation.perturbation;FluidContactEquation.perturbation]));
ContactEquation.perturbation_vector=[ContactEquation.perturbation_vector;FluidContactEquation.perturbation_vector];
ContactEquation.interface_storage=[sparse(maxdof,maxdof),FluidContactEquation.interface_storage.';FluidContactEquation.interface_storage,sparse(maxdof/3,maxdof/3)];
ContactEquation.interface_storage_vector=ContactEquation.interface_storage*[stepU;P];
end