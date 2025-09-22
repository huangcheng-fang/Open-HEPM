%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Efieldcell,CDmatrix]=Perform_constitutive_integration(Cdstrain,Cspin,Efieldcell,CMaterial_id,Material,GNL)
tic;
Cstress=Efieldcell.Stress;
Cstrain=Efieldcell.Strain;
Cpstrain=Efieldcell.Pstrain;
Cparameter=Efieldcell.Parameter;
NumFEcell=numel(Cstress);
%=============================Stress spin==================================
if GNL
for fi=1:1:NumFEcell
    dUdX=Cspin{fi};w=(dUdX-pagetranspose(dUdX))/2;
    R=eye(3)+pagemldivide((eye(3)-0.5*w),w);
    sigmar=reshape(Cstress{fi}([1;4;6;4;2;5;6;5;3],:,:,:),3,3,size(R,3),[]);

    H=Dyadic_R2R2(pageinv(2*eye(3)-w),eye(3)+R,[1,4,2,3]);H=Tensor2matrix_R4(H);
    M=Dyadic_R2R2(eye(3)-dUdX/2,eye(3)+dUdX/2,[1,4,2,3]);M=Tensor2matrix_R4(M);
    Mskew=(M-M([1,4,7,2,5,8,3,6,9],:,:,:))/2;
    S=Dyadic_R2R2(eye(3),pagemtimes(sigmar,'none',R,'transpose'),[1,4,2,3])+...
        Dyadic_R2R2(eye(3),pagemtimes(R,sigmar),[3,1,2,4]);S=Tensor2matrix_R4(S);
    CDrmatrix{fi}=pagemtimes(pagemtimes(S,H),Mskew);
    Ms{fi}=M;
    sigmar=pagemtimes(pagemtimes(R,sigmar),'none',R,'transpose');
    Cstress{fi}=[sigmar(1,1,:,:);sigmar(2,2,:,:);sigmar(3,3,:,:);sigmar(1,2,:,:);sigmar(2,3,:,:);sigmar(1,3,:,:);];
end
end
%======================Constitutive_integration============================
CDmatrix=cell(NumFEcell,1);
for fi=1:1:NumFEcell
    mat=Material(CMaterial_id{fi});
    Estress=Cstress{fi};
    Epstrain=Cpstrain{fi};
    Edstrain=Cdstrain{fi};
    Eparameter=Cparameter{fi};
    
    Cstrain{fi}=Cstrain{fi}+Edstrain;
    [Cstress{fi},Cpstrain{fi},Cparameter{fi},CDmatrix{fi}]=Material_constitutive(mat,Estress,Epstrain,Eparameter,Edstrain);
end
%=============Finite deformation consistent Constitutive matrix============
L=[1,0,0,0,0,0,0,0,0;
   0,0,0,0,1,0,0,0,0;
   0,0,0,0,0,0,0,0,1;
   0,1,0,1,0,0,0,0,0;
   0,0,0,0,0,1,0,1,0;
   0,0,1,0,0,0,1,0,0;];
if GNL
for fi=1:1:NumFEcell
    mat=Material(CMaterial_id{fi});
    Dep=CDmatrix{fi};
    De=Linear_elastic(mat.elasticity.parameter,zeros(6,0),zeros(6,0));
    Dr=CDrmatrix{fi};

    De=pagemrdivide(Dep,De);
    De=De([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3],:,:);
    De([4,8,3],[2,6,7],:,:)=0;   De([2,6,7],[4,8,3],:,:)=0;
    Dr=pagemtimes(De,Dr);

    sigmar=reshape(Cstress{fi}([1;4;6;4;2;5;6;5;3],:,:,:),3,3,size(R,3),[]);
    Ds=Dyadic_R2R2(eye(3),sigmar,[3,4,1,2])-Dyadic_R2R2(eye(3),sigmar,[3,2,1,4]);Ds=Tensor2matrix_R4(Ds);

    Dep=pagemtimes(L.',pagemtimes(pagemtimes(Dep,L),Ms{fi}));

    CDmatrix{fi}=Dep+Dr+Ds;
%     CDmatrix{fi}=pagemtimes(L.',pagemtimes(Dep,L));
end
else
    for fi=1:1:NumFEcell
        Dep=CDmatrix{fi};
        CDmatrix{fi}=pagemtimes(L.',pagemtimes(Dep,L));
    end
end
%==============================Output======================================
Efieldcell.Stress=Cstress;
Efieldcell.Pstrain=Cpstrain;
Efieldcell.Strain=Cstrain;
Efieldcell.Parameter=Cparameter;
%==========================================================================
time=toc;fprintf("Material constitutive integration is completed:%fs\n",time);
end
