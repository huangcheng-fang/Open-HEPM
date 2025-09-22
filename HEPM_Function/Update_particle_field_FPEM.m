%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Pfield,Dmatrix]=Update_particle_field_FPEM(CEdstrain,CEspin,Pfield,Fai,Material,MID,GNL)
tic;
Key=Keywords();
Nstress=Pfield.Stress;
Nstrain=Pfield.Strain;
Npstrain=Pfield.Pstrain;
Nparameter=Pfield.Parameter;
FEcell_num=numel(CEdstrain);
particle_num=size(Nstrain,4);
IP_num=1;
%===========================particle strain================================
Edstrain=[];
for fi=1:FEcell_num
    Esize=size(CEdstrain{fi},[1,2,3,4]);
    Edstrain=[Edstrain;reshape(CEdstrain{fi},Esize(1),[]).'];
end
Ndstrain=reshape((Fai*Edstrain).',6,1,IP_num,particle_num);
%========================particle spin tensor==============================
if GNL
Espin=[];
for fi=1:1:FEcell_num
    Esize=size(CEspin{fi},[1,2,3,4]);
    Espin=[Espin;reshape(CEspin{fi},Esize(1)*Esize(2),[]).'];
end
Nspin=reshape((Fai*Espin)',3,3,IP_num,particle_num);
end
%=============================stress spin==================================
if GNL
    w=eye(3)+pagemldivide((eye(3)-0.5*Nspin),Nspin);
    sigmar=reshape(Nstress([1;4;6;4;2;5;6;5;3],:,:,:),3,3,1,particle_num);
    sigmar=pagemtimes(pagemtimes(w,sigmar),'none',w,'transpose');
    sigmar=reshape(sigmar,9,1,IP_num,particle_num);
    Nstress=sigmar([1;5;9;2;6;3],:,:,:);
end
%=======================Material constitutive==============================
Dmatrix=zeros(6,6,IP_num,particle_num);
for mi=1:1:numel(Material)
nid=find(MID==mi);
mat=Material(mi);
elasticity=mat.elasticity;
plasticity=mat.plasticity;
elasPara=elasticity.parameter;
Mstress=Nstress(:,:,:,nid);
Mpstrain=Npstrain(:,:,:,nid);
Mdstrain=Ndstrain(:,:,:,nid);
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
                [MDep(:,:,ii),Mstress(:,ii),Mpstrain(:,ii)]=Von_Mises_new(elasPara,plasPara,hardPara,[],stress,pstrain,dstrain);
            end
        case Key.mohr_coulomb
            for ii=1:1:numel(nid)*IP_num
                stress=Mstress(:,ii);
                dstrain=Mdstrain(:,ii);
                pstrain=Mpstrain(:,ii);
                [MDep(:,:,ii),Mstress(:,ii),Mpstrain(:,ii)]=Mohr_Coulomb3D_Forward(elasPara,plasPara*pi/180,hardPara,[],stress,pstrain,dstrain);
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
end
%==============================Output======================================
Pfield.Stress=Nstress;
Pfield.Pstrain=Npstrain;
Pfield.Strain=Nstrain+Ndstrain;
Pfield.Parameter=Nparameter;
%==========================================================================
time=toc;fprintf("Element fields is updated:%fs\n",time);
end