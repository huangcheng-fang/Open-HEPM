%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Result,FE]=Update_element_fluid_field(Result,FE,Material,GNL)
tic;
Ev=Result.Edvstrain;
stepU=Result.stepU;
%========================element conductivity matrix=======================
for fi=1:1:numel(FE)
    mat=Material(FE(fi).materialID); 
    FE(fi).D_conductivity=mat.fluid.C_matrix;
    FE(fi).fluid_weight=mat.fluid.weight;
end
%=========================element volume strain============================
if GNL
    for fi=1:1:numel(FE)
        if isempty(FE(fi).D_conductivity);continue;end
        eid=FE(fi).eid;
        SDof=FE(fi).dof;
        dNdX_matrix=FE(fi).F_matrix;%%%%%
        NdU=reshape(stepU(SDof),3,size(SDof,1)/3,1,[]);
        e=-pagemtimes(dNdX_matrix,'none',NdU,'transpose');
        e(1,1,:,:)=e(1,1,:,:)+1;e(2,2,:,:)=e(2,2,:,:)+1;e(3,3,:,:)=e(3,3,:,:)+1;
        Evstrain=-PageDet_mex(e)+1;
        Ev(1,1,1:size(Evstrain,3),eid)=Evstrain;
    end
else
    for fi=1:1:numel(FE)
        if isempty(FE(fi).D_conductivity);continue;end
        eid=FE(fi).eid;
        ndof=size(FE(fi).dof,1);
        EstepU=stepU(reshape(FE(fi).dof,ndof,1,1,[]));
        Evstrain=sum(pagemtimes(FE(fi).B_matrix(1:3,:,:,:),EstepU),1);
        Ev(1,1,1:size(Evstrain,3),eid)=Evstrain;
    end
end
Result.Edvstrain=Ev;
%==========================================================================
time=toc;fprintf("Fluid fields is updated:%fs\n",time);
end