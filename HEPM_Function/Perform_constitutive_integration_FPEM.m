%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Pfield,Dmatrix]=Perform_constitutive_integration_FPEM(dstrain,spin,Pfield,Material,MID,GNL)
tic;
Stress=Pfield.Stress;
Strain=Pfield.Strain;
Pstrain=Pfield.Pstrain;
Parameter=Pfield.Parameter;
particle_num=size(Stress,4);
%=============================Stress spin==================================
if GNL
    w=eye(3)+pagemldivide((eye(3)-0.5*spin),spin);
    sigmar=reshape(Stress([1;4;6;4;2;5;6;5;3],:,:,:),3,3,1,particle_num);
    sigmar=pagemtimes(pagemtimes(w,sigmar),'none',w,'transpose');
    sigmar=reshape(sigmar,9,1,1,particle_num);
    Stress=sigmar([1;5;9;2;6;3],:,:,:);
end
%=======================Material constitutive==============================
Dmatrix=zeros(6,6,1,particle_num);
for mi=1:1:numel(Material)
    nid=find(MID==mi);
    mat=Material(mi);
    Mstress=Stress(:,:,:,nid);
    Mpstrain=Pstrain(:,:,:,nid);
    Mdstrain=dstrain(:,:,:,nid);
    Mparameter=Parameter(:,:,:,nid);
    [Mstress,Mpstrain,Mparameter,MDep]=Material_constitutive(mat,Mstress,Mpstrain,Mparameter,Mdstrain);
    Dmatrix(:,:,:,nid)=MDep;
    Stress(:,:,:,nid)=Mstress;
    Pstrain(:,:,:,nid)=Mpstrain;
    Parameter(:,:,:,nid)=Mparameter;
end
%==============================Output======================================
Pfield.Stress=Stress;
Pfield.Strain=Strain+dstrain;
Pfield.Pstrain=Pstrain;
Pfield.Parameter=Parameter;
%==========================================================================
time=toc;fprintf("Material constitutive integration is completed:%fs\n",time);
end