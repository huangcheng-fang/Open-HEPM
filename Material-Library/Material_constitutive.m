%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Estress,Epstrain,Eparameter,EDep]=Material_constitutive(mat,Estress,Epstrain,Eparameter,Edstrain)
Key=Keywords();
elasticity=mat.elasticity;
plasticity=mat.plasticity;
elasPara=elasticity.parameter;
stress_num=size(Estress,1);
ip_num=size(Estress,3);
element_num=size(Estress,4);
%-------------------------elastic constitutive-----------------------------
if plasticity.type==Key.none
    switch elasticity.type
        case Key.linear
            [EDep,Estress]=Linear_elastic(elasPara,Edstrain,Estress);
            EDep=repmat(EDep,1,1,size(Edstrain,3),size(Edstrain,4));
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
    EDep=zeros(stress_num,stress_num,ip_num,element_num);
    switch plasticity.type
        case Key.drucker_prager
            for ii=1:1:element_num
                for ipi=1:1:ip_num
                    stress=Estress(:,1,ipi,ii);
                    dstrain=Edstrain(:,1,ipi,ii);
                    pstrain=Epstrain(:,1,ipi,ii);
                    [EDep(:,:,ipi,ii),Estress(:,1,ipi,ii),Epstrain(:,1,ipi,ii)]=Easy_Drucker_Prager(elasticity.parameter,plasPara,hardPara,[],stress,pstrain,dstrain);
                end
            end
        case Key.mohr_coulomb
            for ii=1:1:element_num
                for ipi=1:1:ip_num
                    stress=Estress(:,1,ipi,ii);
                    dstrain=Edstrain(:,1,ipi,ii);
                    pstrain=Epstrain(:,1,ipi,ii);
                    [EDep(:,:,ipi,ii),Estress(:,1,ipi,ii),Epstrain(:,1,ipi,ii)]=Mohr_Coulomb3D(elasPara,plasPara/180*pi,hardPara,[],stress,pstrain,dstrain);
                end
            end
        case Key.tresca
            for ii=1:1:element_num
                for ipi=1:1:ip_num
                    stress=Estress(:,1,ipi,ii);
                    dstrain=Edstrain(:,1,ipi,ii);
                    pstrain=Epstrain(:,1,ipi,ii);
                    [Estress(:,1,ipi,ii),EDep(:,:,ipi,ii),Epstrain(:,1,ipi,ii)]=Easy_Tresca(elasticity.parameter,plasPara,[],hardPara,stress,pstrain,dstrain);
                end
            end
        case Key.mises
            for ii=1:1:element_num
                for ipi=1:1:ip_num
                    stress=Estress(:,1,ipi,ii);
                    dstrain=Edstrain(:,1,ipi,ii);
                    pstrain=Epstrain(:,1,ipi,ii);
                    histPara=Eparameter(1:2,1,ipi,ii);
                    % [EDep(:,:,ipi,ii),Estress(:,1,ipi,ii),Epstrain(:,1,ipi,ii)]=Von_Mises_new(elasPara,plasPara,hardPara,[],stress,pstrain,dstrain);
                    [Estress(:,1,ipi,ii),EDep(:,:,ipi,ii),Epstrain(:,1,ipi,ii),Eparameter(1:2,1,ipi,ii)]=Easy_Mises(elasticity.parameter,plasPara,histPara,[],stress,pstrain,dstrain);
                end
            end
        case Key.modified_cam_clay
            if elasticity.type~=Key.porous_elastic;error('Modified_cam_clay must use porous elastic material');end
            if size(Eparameter,1)<2;error('Please use Creat_predefined_field to define initial void ratio and consolidation pressure');end
            error('MCC')
        otherwise
            error('Unknown plastic material')
    end
end
end