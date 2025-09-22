%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Contact_pair=Get_contact_stress_new(Contact_pair,Contact_condition,Fcon,positions,stepU)
coor=positions';coor=coor(:);
%=========================contact stress===================================
for ci=1:1:numel(Contact_pair)
    R_matrix=Contact_pair(ci).R_matrix;
    contact_matrix=Contact_pair(ci).contact_matrix;
    contact_status=Contact_pair(ci).contact_status;
    slave_dof=Contact_pair(ci).slave_dof;
    area_vector=Contact_pair(ci).area_vector;
    contact_stress=R_matrix*Contact_pair(ci).NCstress(slave_dof)./area_vector;
    Em=Contact_condition(ci).parameter.Emodulus;
    switch Contact_condition(ci).enforcement
        case 'dual_lagrange'
            contact_stress=R_matrix*(Fcon(slave_dof)./area_vector);
        case 'augmented_lagrange'
            % dcontact_stress=Em*(R_matrix*(contact_matrix*stepU));
            % dcontact_stress(3:3:end)=Em*(R_matrix(3:3:end,:)*(contact_matrix*coor));
            % contact_stress(contact_status)=contact_stress(contact_status)+dcontact_stress(contact_status);
            % contact_stress(contact_status)=contact_stress(contact_status)+Em*(R_matrix(contact_status,:)*(contact_matrix*coor));
            contact_stress(contact_status)=R_matrix(contact_status,:)*(Fcon(slave_dof)./area_vector);
%             contact_stress=R_matrix*(Fcon(slave_dof)./area_vector);
        otherwise
            error(['Unknow contact enforcement: ',Contact_condition(ci).enforcement])
    end
    Contact_pair(ci).contact_stress=contact_stress;
end
end