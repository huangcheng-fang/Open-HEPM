%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function ContactEquation=Update_contact_Equation_new(ContactEquation,Contact_pair,Contact_condition,maxdof)
equation_change=false;
for ci=1:1:numel(Contact_pair)
    if Contact_pair(ci).update||Contact_pair(ci).status_change
        equation_change=true;
    end
end
ContactEquation.update=equation_change;
%==========================================================================
Penalty.matrix=sparse(maxdof,maxdof);
Penalty.vector=zeros(maxdof,1);
Dual_Lagrange.matrix=speye(maxdof);
Dual_Lagrange.Fcon=zeros(maxdof,1);
Dual_Lagrange.slave_dof=zeros(0,1);
Dual_Lagrange.constraint=zeros(0,2);
Dual_Lagrange.perturbation=sparse(maxdof,maxdof);
for ci=1:1:numel(Contact_pair)
    R_matrix=Contact_pair(ci).R_matrix;
    slave_dof=Contact_pair(ci).slave_dof;
    contact_gap=Contact_pair(ci).contact_gap;
    gapT1=contact_gap(1:3:end);
    gapT2=contact_gap(2:3:end);
    area_vector=Contact_pair(ci).area_vector;
    contact_matrix=Contact_pair(ci).contact_matrix;
    contact_stress=Contact_pair(ci).contact_stress;
    contact_status=Contact_pair(ci).contact_status;
    contact_statusT=contact_status(1:3:end);
    contact_statusN=contact_status(3:3:end);
    mu=Contact_condition(ci).parameter.friction_coeff;
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
            
            perturbation=(area_vector(3:3:end).*[Gm,Gm,Em]).';
            perturbation(~contact_status)=0;
            perturbation=sparse(slave_dof,slave_dof,perturbation,maxdof,maxdof);
            Dual_Lagrange.perturbation=Dual_Lagrange.perturbation+perturbation;
        case {'augmented_lagrange'}
            RContact=R_matrix(contact_status,:)*contact_matrix;
            contact_status0=contact_status;
            % contact_status0(1:3:end)=0;contact_status0(2:3:end)=0;
            contact_stress(contact_status0)=Em*contact_gap(contact_status0);
            Penalty.vector=Penalty.vector+contact_matrix.'*(R_matrix.'*(area_vector.*contact_stress));
            Penalty.matrix=Penalty.matrix+RContact.'*(Em.*area_vector(contact_status).*RContact);
        otherwise
            error(['Unknow contact enforcement: ',Contact_condition(ci).enforcement])
    end
end
%============================Output========================================
ContactEquation.Penalty=Penalty;
ContactEquation.Dual_Lagrange=Dual_Lagrange;
end