%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Contact_pair=Get_contact_stress(Contact_pair,Contact_condition,Fcon,Lamda)
flag=1;
for ci=1:1:numel(Contact_pair)
    R_matrix=Contact_pair(ci).R_matrix;
    %contact_status=Contact_pair(ci).contact_status;
    slave_dof=Contact_pair(ci).slave_dof;
    area_vector=Contact_pair(ci).area_vector;
    switch Contact_condition(ci).enforcement
        case 'dual_lagrange'
            contact_stress=R_matrix*(Fcon(slave_dof)./area_vector);
            Contact_pair(ci).contact_stress_vector(slave_dof,:)=contact_stress;
        case 'lagrange'
            contact_stress=Lamda(flag:flag+numel(slave_dof)-1);
            Contact_pair(ci).contact_stress_vector(slave_dof,:)=Contact_pair(ci).contact_stress_vector(slave_dof,:)+contact_stress;
            flag=flag+numel(slave_dof);
        otherwise
            error(['Unknow contact enforcement: ',Contact_pair(ci).enforcement])
    end
end
%==========================================================================
for ci=1:1:numel(Contact_pair)
    %R_matrix=Contact_pair(ci).R_matrix;
    %contact_status=Contact_pair(ci).contact_status;
    slave_nid=Contact_pair(ci).slave_dof(3:3:end)/3;
    %area_vector=Contact_pair(ci).area_vector;
    switch Contact_condition(ci).enforcement
        case 'dual_lagrange'
            error('')
        case 'lagrange'
            contact_flux=Lamda(flag:flag+numel(slave_nid)-1);
            Contact_pair(ci).contact_flux_vector(slave_nid,:)=Contact_pair(ci).contact_flux_vector(slave_nid,:)+contact_flux;
            flag=flag+numel(slave_nid);
        otherwise
            error(['Unknow contact enforcement: ',Contact_pair(ci).enforcement])
    end
end
end