%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Pfield,Dmatrix]=Update_particle_field(CEdstrain,CEspin,Pfield,Fai,FEcell,Material,MID,GNL)
tic;
Key=Keywords();
Nstress=Pfield.Stress;
Nstrain=Pfield.Strain;
Npstrain=Pfield.Pstrain;
Nparameter=Pfield.Parameter;
CR_matrix=FEcell.R_matrix;
CEid=FEcell.Element_id;
FEcell_num=numel(CR_matrix);
node_num=size(Nstrain,4);
element_num=size(Fai{1},2);
IP_num=numel(Fai);
%============================element strain================================

%===========================particle strain================================
Ndstrain=zeros(size(Nstrain,1),1,IP_num,node_num);
for fi=1:1:FEcell_num
    Edstrain=zeros(element_num,size(CEdstrain{fi},1),1,size(CEdstrain{fi},3));
    Edstrain(CEid{fi},:,:,:)=permute(CEdstrain{fi},[4,1,2,3]);
    for ipi=1:1:IP_num
        Ndstrain(:,1,ipi,:)=Ndstrain(:,1,ipi,:)+permute(Fai{ipi}*Edstrain,[2,3,4,1]);
    end
end
%========================particle spin tensor==============================
Nspin=zeros(9,1,IP_num,node_num);
if GNL
for fi=1:1:FEcell_num
    Espin=zeros(element_num,size(CEspin{fi},1)*size(CEspin{fi},2),1,size(CEspin{fi},3));
    Espin(CEid{fi},:,:,:)=reshape(permute(CEspin{fi},[4,1,2,3]),[],size(Espin,2),1,size(Espin,4));
    for ipi=1:1:IP_num
        Nspin(:,1,ipi,:)=Nspin(:,1,ipi,:)+permute(Fai{ipi}*Espin,[2,3,4,1]);
    end
end
end
%=============================stress spin==================================
if GNL
for ipi=1:1:IP_num
    w=reshape(Nspin(:,:,ipi,:),3,3,1,node_num);
    w=eye(3)+pagemldivide((eye(3)-0.5*w),w);
    sigmar=reshape(Nstress([1;4;6;4;2;5;6;5;3],:,ipi,:),3,3,1,node_num);
    sigmar=pagemtimes(pagemtimes(w,sigmar),'none',w,'transpose');
    Nstress(:,:,ipi,:)=[sigmar(1,1,:,:);sigmar(2,2,:,:);sigmar(3,3,:,:);sigmar(1,2,:,:);sigmar(2,3,:,:);sigmar(1,3,:,:);];
end
end
%=======================Material constitutive==============================
Dmatrix=zeros(6,6,IP_num,node_num);
for mi=1:1:numel(Material)
nid=find(MID==mi);
mat=Material(mi);
elasticity=mat.elasticity;
plasticity=mat.plasticity;
elasPara=elasticity.parameter;
Mstress=Nstress(:,:,:,nid);
Mpstrain=Npstrain(:,:,:,nid);
Mdstrain=Ndstrain(:,:,:,nid);
Mparameter=Nparameter(:,:,:,nid);
%-------------------------elastic constitutive-----------------------------
if plasticity.type==Key.none
    switch elasticity.type
        case Key.linear
            [MDep,Mstress]=Linear_elastic(elasPara,Mdstrain,Mstress);
            MDep=repmat(MDep,1,1,IP_num,numel(nid));
        otherwise
            error('Unknown elastic material')
    end
end
%--------------------------plastic constitutive----------------------------
if plasticity.type~=Key.none
    if elasticity.type==Key.linear
        elasPara=Linear_elastic(elasPara,zeros(6,0),zeros(6,0));
    end
    plasPara=plasticity.parameter;
    hardPara=plasticity.harden_parameter;
    MDep=zeros(6,6,IP_num,numel(nid));
    switch plasticity.type
        case Key.mises
            for ii=1:1:numel(nid)*IP_num
                stress=Mstress(:,ii);
                dstrain=Mdstrain(:,ii);
                pstrain=Mpstrain(:,ii);
                hisPara=Mparameter(:,ii);
%                 [MDep(:,:,ii),Mstress(:,ii),Mpstrain(:,ii)]=Von_Mises_new(elasPara,plasPara,hardPara,[],stress,pstrain,dstrain);
                [Mstress(:,ii),MDep(:,:,ii),Mpstrain(:,ii),Mparameter(:,ii)]=Easy_Mises(elasticity.parameter,plasPara,hisPara,hardPara,stress,pstrain,dstrain);
            end
        case Key.mohr_coulomb
            for ii=1:1:numel(nid)*IP_num
                stress=Mstress(:,ii);
                dstrain=Mdstrain(:,ii);
                pstrain=Mpstrain(:,ii);
                [MDep(:,:,ii),Mstress(:,ii),Mpstrain(:,ii)]=Mohr_Coulomb3D_Forward(elasPara,plasPara*pi/180,hardPara,[],stress,pstrain,dstrain);
            end
        case Key.tresca
            for ii=1:1:numel(nid)*IP_num
                stress=Mstress(:,ii);
                dstrain=Mdstrain(:,ii);
                pstrain=Mpstrain(:,ii);
                [Mstress(:,ii),MDep(:,:,ii),Mpstrain(:,ii)]=Easy_Tresca(elasticity.parameter,plasPara,[],hardPara,stress,pstrain,dstrain);
            end
        case Key.drucker_prager
            for ii=1:1:numel(nid)*IP_num
                stress=Mstress(:,ii);
                dstrain=Mdstrain(:,ii);
                pstrain=Mpstrain(:,ii);
                [MDep(:,:,ii),Mstress(:,ii),Mpstrain(:,ii)]=Drucker_Prager3D(elasPara,plasPara,hardPara,[],stress,pstrain,dstrain);
            end
        otherwise
            error('Unknown plastic material')
    end
end
%--------------------------------------------------------------------------
Dmatrix(:,:,:,nid)=MDep;
Nstress(:,:,:,nid)=Mstress;
Npstrain(:,:,:,nid)=Mpstrain;
Nparameter(:,:,:,nid)=Mparameter;
end
%==============================Output======================================
Pfield.Stress=Nstress;
Pfield.Pstrain=Npstrain;
Pfield.Strain=Nstrain+Ndstrain;
Pfield.Parameter=Nparameter;
%==========================================================================
time=toc;fprintf("Element fields is updated:%fs\n",time);
end