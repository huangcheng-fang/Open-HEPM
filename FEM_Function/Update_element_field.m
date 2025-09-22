%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Efieldcell,CDmatrix]=Update_element_field(Cdstrain,Cspin,Efieldcell,FEcell,Material,GNL)
tic;
Key=Keywords();
Cstress=Efieldcell.Stress;
Cstrain=Efieldcell.Strain;
Cpstrain=Efieldcell.Pstrain;
Cparameter=Efieldcell.Parameter;
CR_matrix=FEcell.R_matrix;
CEid=FEcell.Element_id;
CMaterial_id=FEcell.Material_id;
NumFEcell=numel(CR_matrix);
%=============================stress spin==================================
if GNL
    for fi=1:1:NumFEcell
        if ~isempty(CR_matrix{fi}); continue; end
        w=Cspin{fi};
        w=eye(3)+pagemldivide((eye(3)-0.5*w),w);
        sigmar=reshape(Cstress{fi}([1;4;6;4;2;5;6;5;3],:,:,:),3,3,size(w,3),[]);
        sigmar=pagemtimes(pagemtimes(w,sigmar),'none',w,'transpose');
        Cstress{fi}=[sigmar(1,1,:,:);sigmar(2,2,:,:);sigmar(3,3,:,:);sigmar(1,2,:,:);sigmar(2,3,:,:);sigmar(1,3,:,:);];
    end
end
%=======================Material constitutive==============================
CDmatrix=cell(NumFEcell,1);
for fi=1:1:NumFEcell
eid=CEid{fi};
mat=Material(CMaterial_id{fi});
elasticity=mat.elasticity;
plasticity=mat.plasticity;
elasPara=elasticity.parameter;
Estress=Cstress{fi};
Epstrain=Cpstrain{fi};
Edstrain=Cdstrain{fi};
stress_num=size(Estress,1);
ip_num=size(Estress,3);
element_num=numel(eid);
%-------------------------elastic constitutive-----------------------------
if plasticity.type==Key.none
    switch elasticity.type
        case Key.linear
            [EDep,Estress]=Linear_elastic(elasPara,Edstrain,Estress);
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
        case Key.mohr_coulomb
            for ii=1:1:element_num
                for ipi=1:1:ip_num
                    stress=Estress(:,1,ipi,ii);
                    dstrain=Edstrain(:,1,ipi,ii);
                    pstrain=Epstrain(:,1,ipi,ii);
                    [EDep(:,:,ipi,ii),Estress(:,1,ipi,ii),Epstrain(:,1,ipi,ii)]=Mohr_Coulomb3D(elasPara,plasPara/180*pi,hardPara,[],stress,pstrain,dstrain);
                end
            end
        case Key.mises
            for ii=1:1:element_num
                for ipi=1:1:ip_num
                    stress=Estress(:,1,ipi,ii);
                    dstrain=Edstrain(:,1,ipi,ii);
                    pstrain=Epstrain(:,1,ipi,ii);
                    [EDep(:,:,ipi,ii),Estress(:,1,ipi,ii),Epstrain(:,1,ipi,ii)]=Von_Mises_new(elasPara,plasPara,hardPara,[],stress,pstrain,dstrain);
                end
            end
        case Key.modified_cam_clay
            if elasticity.type~=Key.porous_elastic;error('Modified_cam_clay must use porous elastic material');end
            if size(Cparameter,1)<2;error('Please use Creat_predefined_field to define initial void ratio and consolidation pressure');end
            for ii=1:1:numel(eid)
                ii=eid(ii);
                for ipi=1:1:ip_num
                    Estress=Cstress(:,1,ipi,ii);
                    Edstrain=Cdstrain(:,1,ipi,ii);
                    pstrain=Cpstrain(:,1,ipi,ii);
                    voidratio=Cparameter(1,1,ipi,ii);
                    pc=Cparameter(2,1,ipi,ii);
                    [Cstress(:,1,ipi,ii),Cpstrain(:,1,ipi,ii),Cparameter(1:2,1,ipi,ii),EDep(:,:,ipi,ii)]=Modified_Cam_Clay3D_new(elasPara,plasPara,[],Estress,pstrain,Edstrain,[voidratio,pc]);
                end
            end
        otherwise
            error('Unknown plastic material')
    end
end
%--------------------------------------------------------------------------
CDmatrix{fi}=EDep;
Cstress{fi}=Estress;
Cstrain{fi}=Cstrain{fi}+Edstrain;
Cpstrain{fi}=Epstrain;
end
%==============================Output======================================
Efieldcell.Stress=Cstress;
Efieldcell.Pstrain=Cpstrain;
Efieldcell.Strain=Cstrain;
Efieldcell.Parameter=Cparameter;
%==========================================================================
time=toc;fprintf("Element fields is updated:%fs\n",time);
end