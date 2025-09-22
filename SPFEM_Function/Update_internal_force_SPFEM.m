%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Fint=Update_internal_force_SPFEM(Pstress,Fai,node_volume,FEcell,FEUdof,maxdof)
tic;
CEid=FEcell.Element_id;
CDetJac=FEcell.DetJac;
CB_matrix=FEcell.B_matrix;
NumFEcell=numel(CDetJac);
%==========================================================================
volume_matrix=diag(node_volume/6);
Estress=0;
for ipi=1:1:numel(Fai)
Estress=Estress+permute(Pstress(:,1,ipi,:),[1,4,3,2])*(volume_matrix*Fai{ipi});
end
Estress=permute(Estress,[1,4,3,2]);
%==========================================================================
Fint=zeros(maxdof,1);
for fi=1:1:NumFEcell
    stress=Estress(:,:,:,CEid{fi});
    dof=FEUdof{fi};
    Fint_element=pagemtimes(CB_matrix{fi},'transpose',stress,'none');
    % Fint_element=sum(pagemtimes(Fint_element,CDetJac{fi}),3);
    Fint=Fint+accumarray(dof(:),Fint_element(:),[maxdof,1]);
end
%==========================================================================
time=toc;fprintf("Internal force is updated:%fs\n",time);
end