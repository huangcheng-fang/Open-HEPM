%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Result]=Get_contact_stress(Result,Contact_pair,Contact_condition,Fcon,Lamda)
dU=Result.dU;
NCstress=Result.NCstress;
NCflux=Result.NCflux;
maxdof=size(dU,1);
flag=1;
%=========================contact stress===================================
for ci=1:1:numel(Contact_pair)
    R_matrix=Contact_pair(ci).R_matrix;
    slave_dof=Contact_pair(ci).slave_dof;
    area_vector=Contact_pair(ci).area_vector;
    contact_matrix=Contact_pair(ci).contact_matrix;
    contact_stress=Contact_pair(ci).contact_stress;
    contact_status=Contact_pair(ci).contact_status;
    En=Contact_condition(ci).parameter.Emodulus;
    switch Contact_condition(ci).enforcement
        case 'dual_lagrange'
            contact_stress=R_matrix*(Fcon(slave_dof)./area_vector);
            NCstress(slave_dof,ci)=contact_stress;
        case {'perturbed_lagrange','lagrange'}
            NCstress(slave_dof,ci)=contact_stress;
            active_dof=slave_dof(contact_status);
            contact_stress=Lamda(flag:flag+numel(active_dof)-1);
            NCstress(active_dof,ci)=contact_stress+NCstress(active_dof,ci);
            flag=flag+numel(active_dof);
        case 'penalty'
            NCstress(slave_dof,ci)=contact_stress;
            active_dof=slave_dof(contact_status);
            contact_stress=En*R_matrix(contact_status,:)*(contact_matrix*dU);
            NCstress(active_dof,ci)=contact_stress+NCstress(active_dof,ci);
        otherwise
            error(['Unknow contact enforcement: ',Contact_condition(ci).enforcement])
    end
end
%=====================stabilization term===================================
% inactive_dofs=zeros(size(Fcon,1),1);
% for ci=1:1:numel(Contact_pair)
%     slave_dof=Contact_pair(ci).slave_dof;
%     contact_status=Contact_pair(ci).contact_status;
%     inactive_dofs(slave_dof(~contact_status))=1;   
% end
% stable_multiplier(flag:flag+sum(inactive_dofs,'all')-1)=Lamda(flag:flag+sum(inactive_dofs,'all')-1);
% flag=flag+sum(inactive_dofs,'all');
%==========================contact flux====================================
try
for ci=1:1:numel(Contact_pair)
    slave_nid=Contact_pair(ci).slave_dof(3:3:end)/3;
    area_vector=Contact_pair(ci).area_vector(3:3:end);
    switch Contact_condition(ci).enforcement
        case 'dual_lagrange'
            contact_flux=Fcon(slave_nid+maxdof)./area_vector;
            NCflux(slave_nid,ci)=contact_flux;
        case {'perturbed_lagrange','lagrange'}
            contact_flux=Lamda(flag:flag+numel(slave_nid)-1);
            NCflux(slave_nid,ci)=NCflux(slave_nid,ci)+contact_flux;
            flag=flag+numel(slave_nid);
        case 'penalty'

        otherwise
            error(['Unknow contact enforcement: ',Contact_condition(ci).enforcement])
    end
end
catch
disp('nnnn')
end
%===============================Output=====================================
Result.NCstress=NCstress;
Result.NCflux=NCflux;
end