%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function ContactEquation=Update_contact_Equation(Contact_pair,Contact_condition,maxdof)
ContactEquation.Penalty.matrix=sparse(maxdof,maxdof);
ContactEquation.Penalty.vector=zeros(maxdof,1);
Dual_Lagrange.matrix=speye(maxdof);
Dual_Lagrange.Fcon=zeros(maxdof,1);
Dual_Lagrange.slave_dof=zeros(0,1);
Dual_Lagrange.constraint=zeros(0,2);
Dual_Lagrange.perturbation=sparse(maxdof,1);


ContactEquation.Lagrange=sparse(0,maxdof);
ContactEquation.perturbation=sparse(0,1);
ContactEquation.perturbation_vector=sparse(0,1);
ContactEquation.Fcon=zeros(maxdof,1);
ContactEquation.constraint=zeros(0,2);
stabilization=sparse(maxdof,maxdof);
for ci=1:1:numel(Contact_pair)
    R_matrix=Contact_pair(ci).R_matrix;
    slave_dof=Contact_pair(ci).slave_dof;
    contact_gap=Contact_pair(ci).contact_gap;
    area_vector=Contact_pair(ci).area_vector;
    contact_matrix=Contact_pair(ci).contact_matrix;
    contact_stress=Contact_pair(ci).contact_stress;
    contact_status=Contact_pair(ci).contact_status;
    plane_orientation=Contact_condition(ci).plane_orientation;
    Em=Contact_condition(ci).parameter.Emodulus;
    Gm=Contact_condition(ci).parameter.Gmodulus;
    if strcmp(Contact_condition(ci).constraint,'tie');contact_gap=contact_gap*0;end
    switch Contact_condition(ci).enforcement
        case 'dual_lagrange'
            Dual_Lagrange.Fcon=Dual_Lagrange.Fcon+contact_matrix.'*(R_matrix.'*(area_vector.*contact_stress));
            contact_matrix=-contact_matrix;
            contact_matrix(:,slave_dof)=R_matrix.';
            Dual_Lagrange.matrix(slave_dof,:)=contact_matrix;
            Dual_Lagrange.slave_dof=[Dual_Lagrange.slave_dof;slave_dof];
            if Contact_condition(ci).perturbation
                perturbation=(area_vector(3:3:end).*[Gm,Gm,Em]).';perturbation=perturbation(:);
                Dual_Lagrange.perturbation(slave_dof(contact_status),1)=perturbation(contact_status);
            else
                Dual_Lagrange.constraint=[Dual_Lagrange.constraint;slave_dof(contact_status),-contact_gap(contact_status)];
            end
        case 'lagrange'
            contact_matrix=area_vector.*contact_matrix;
            ContactEquation.Fcon=ContactEquation.Fcon+contact_matrix.'*(R_matrix.'*(contact_stress));
            ContactEquation.Lagrange=[ContactEquation.Lagrange;R_matrix(contact_status,:)*contact_matrix];
            perturbation=zeros(numel(slave_dof),1);
            ContactEquation.perturbation=[ContactEquation.perturbation;-perturbation(contact_status)];
            perturbation_vector=[contact_gap*0,contact_gap*0,contact_gap].';perturbation_vector=perturbation_vector(:).*area_vector;
            ContactEquation.perturbation_vector=[ContactEquation.perturbation_vector;-max(perturbation_vector(contact_status),0)];
            inactive_dof=slave_dof(~contact_status);
%             stabilization(inactive_dof,:)=stabilization(inactive_dof,:)+plane_orientation*R_matrix(~contact_status,:)*contact_matrix;
        case 'perturbed_lagrange'
            contact_matrix=area_vector.*contact_matrix;
            ContactEquation.Fcon=ContactEquation.Fcon+contact_matrix.'*(R_matrix.'*(contact_stress));
            ContactEquation.Lagrange=[ContactEquation.Lagrange;R_matrix(contact_status,:)*contact_matrix];
            perturbation=(max(-contact_gap,1e-8).*area_vector(3:3:end).*[1/Gm,1/Gm,1/Em]).';perturbation=perturbation(:);
            ContactEquation.perturbation=[ContactEquation.perturbation;-perturbation(contact_status)];
            perturbation_vector=[contact_gap*0,contact_gap*0,contact_gap].';perturbation_vector=perturbation_vector(:).*area_vector;
            ContactEquation.perturbation_vector=[ContactEquation.perturbation_vector;-max(perturbation_vector(contact_status),0)];
            inactive_dof=slave_dof(~contact_status);
            stabilization(inactive_dof,:)=stabilization(inactive_dof,:)+plane_orientation*R_matrix(~contact_status,:)*contact_matrix;
        case 'penalty'
            contact_matrix=area_vector.*contact_matrix;
            ContactEquation.Penalty.vector=ContactEquation.Penalty.vector+contact_matrix.'*(R_matrix.'*(contact_stress));
            RContact=R_matrix(contact_status,:)*contact_matrix;
            ContactEquation.Penalty.matrix=ContactEquation.Penalty.matrix+RContact.'*(Em./area_vector(contact_status).*RContact);
%             RContact=R_matrix(~contact_status,:)*contact_matrix;
%             ContactEquation.Penalty.matrix=ContactEquation.Penalty.matrix+RContact.'*(Em/1e8./area_vector(~contact_status).*RContact);
        otherwise
            error(['Unknow contact enforcement: ',Contact_condition(ci).enforcement])
    end
end
%============================Output========================================
ContactEquation.Dual_Lagrange=Dual_Lagrange;
end