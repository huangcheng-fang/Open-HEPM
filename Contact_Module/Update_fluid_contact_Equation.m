%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FluidContactEquation=Update_fluid_contact_Equation(Contact_pair,Contact_condition,P,dU,maxdof)
%==============================Interface_thickness=========================
thickness=zeros(maxdof,1);
for ci=1:numel(Contact_pair)
    slave_nid=Contact_pair(ci).slave_dof(3:3:end)/3;
    thickness(slave_nid)=thickness(slave_nid)+Contact_pair(ci).contact_gap(3:3:end);
end
thickness=max(-thickness,2e-5);
%==========================================================================
FluidContactEquation.Penalty.matrix=sparse(maxdof,maxdof);
FluidContactEquation.Penalty.vector=zeros(maxdof,1);
Dual_Lagrange.matrix=speye(maxdof);
Dual_Lagrange.Fcon=zeros(maxdof,1);
Dual_Lagrange.slave_dof=zeros(0,1);
Dual_Lagrange.constraint=zeros(0,2);
Dual_Lagrange.perturbation=sparse(maxdof,1);

FluidContactEquation.Fcon=zeros(maxdof,1);
FluidContactEquation.Lagrange=sparse(0,maxdof);
FluidContactEquation.perturbation=sparse(0,1);
FluidContactEquation.perturbation_vector=sparse(0,1);
FluidContactEquation.constraint=zeros(0,2);
FluidContactEquation.interface_storage=sparse(maxdof,3*maxdof);
for ci=1:1:numel(Contact_pair)
    slave_dof=Contact_pair(ci).slave_dof(3:3:end)/3;
    area_vector=Contact_pair(ci).area_vector(3:3:end);
    contact_gap=Contact_pair(ci).contact_gap;
    contact_gap_P=Contact_pair(ci).contact_gap_P;
    contact_flux=Contact_pair(ci).contact_flux;
    contact_matrix=Contact_pair(ci).contact_matrix(3:3:end,3:3:end);
    R_contact_matrix=Contact_pair(ci).R_matrix(3:3:end,:)*Contact_pair(ci).contact_matrix;
    permeability=Contact_condition(ci).parameter.permeability;
    switch Contact_condition(ci).enforcement
        case 'dual_lagrange'
            Dual_Lagrange.Fcon=Dual_Lagrange.Fcon+contact_matrix.'*((area_vector.*contact_flux));
            contact_matrix=-contact_matrix;
            contact_matrix(:,slave_dof)=speye(numel(slave_dof));
            Dual_Lagrange.matrix(slave_dof,:)=contact_matrix;
            Dual_Lagrange.slave_dof=[Dual_Lagrange.slave_dof;slave_dof+maxdof*3];
            if Contact_condition(ci).perturbation
                perturbation=area_vector*permeability;
                Dual_Lagrange.perturbation(slave_dof,1)=-perturbation;
            else
                Dual_Lagrange.constraint=[Dual_Lagrange.constraint;slave_dof+maxdof*3,zeros(numel(slave_dof),1)];
            end
        case {'perturbed_lagrange','lagrange'}
            FluidContactEquation.Fcon=FluidContactEquation.Fcon+contact_matrix.'*(area_vector.*contact_flux);
            FluidContactEquation.Lagrange=[FluidContactEquation.Lagrange;area_vector.*contact_matrix];
            storage=area_vector.*R_contact_matrix;
            FluidContactEquation.interface_storage(slave_dof,:)=FluidContactEquation.interface_storage(slave_dof,:)+storage;
            contact_gap(abs(contact_gap)<1e-8)=0;
            perturbation=max(-contact_gap,1e-8).*area_vector*permeability^-1;
            FluidContactEquation.perturbation=[FluidContactEquation.perturbation;perturbation];
            FluidContactEquation.perturbation_vector=[FluidContactEquation.perturbation_vector;-area_vector.*contact_gap_P-perturbation.*contact_flux];
        case 'penalty'
            contact_matrix=area_vector.*contact_matrix;
            contact_thickness=thickness(slave_dof)/2;
            cn=permeability./max(contact_thickness,1e-5)./area_vector;
            FluidContactEquation.Penalty.matrix=FluidContactEquation.Penalty.matrix-contact_matrix.'*(cn.*contact_matrix);
            FluidContactEquation.Penalty.vector=FluidContactEquation.Penalty.vector-contact_matrix.'*(cn.*contact_matrix*P);
            storage=area_vector.*R_contact_matrix;
            FluidContactEquation.interface_storage(slave_dof,:)=FluidContactEquation.interface_storage(slave_dof,:)+storage;
        otherwise
            error(['Unknow contact enforcement: ',Contact_condition(ci).enforcement])
    end
end
%============================Output========================================
FluidContactEquation.Dual_Lagrange=Dual_Lagrange;
end