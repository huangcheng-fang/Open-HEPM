%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Fint=Update_internal_force_FPEM(Pstress,Fai,node_volume,FEcell,FEUdof,maxdof)
tic;
CEid=FEcell.Element_id;
CB_matrix=FEcell.B_matrix;
NumFEcell=numel(CEid);
%==========================================================================
volume_matrix=diag(node_volume);
Estress=permute(Pstress,[1,4,3,2])*(volume_matrix*Fai);
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