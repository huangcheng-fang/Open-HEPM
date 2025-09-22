%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [CPdstrain,CPspin]=Get_particle_dstrain(CEdstrain,CEspin,CEMid,FPcell,GNL)
tic;
CSmoothing_matrix=FPcell.Smoothing_matrix;
CPMid=FPcell.Material_id;
FPcell_num=numel(CSmoothing_matrix);
%=============================Collect result===============================
CEdstrain=Collect_cell4D(CEdstrain,CEMid,CPMid);
CEspin=Collect_cell4D(CEspin,CEMid,CPMid);
%=========================particle strain and spin tensor==================
CPdstrain=coder.nullcopy(cell(FPcell_num,1));
CPspin=coder.nullcopy(cell(FPcell_num,1));
for pci=1:1:FPcell_num
    Smoothing_matrix=CSmoothing_matrix{pci};
    particle_num=size(Smoothing_matrix,1);
    Esize=size(CEdstrain{pci},[1,2,3,4]);
    Edstrain=reshape(CEdstrain{pci},Esize(1)*Esize(2),[]).';
    CPdstrain{pci}=reshape((Smoothing_matrix*Edstrain).',Esize(1),Esize(2),1,particle_num);
    %========================particle spin tensor==============================
    if ~GNL; continue; end
    Esize=size(CEspin{pci},[1,2,3,4]);
    Espin=reshape(CEspin{pci},Esize(1)*Esize(2),[]).';
    CPspin{pci}=reshape((Smoothing_matrix*Espin)',Esize(1),Esize(2),1,particle_num);
end
%==========================================================================
time=toc;fprintf("Particle dstrain is updated:%fs\n",time);
end