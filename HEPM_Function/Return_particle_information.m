%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [particles,Pinpart]=Return_particle_information(FPcell,Part_list)
Particles=FPcell.Particles;
%==========================================================================
particles=Montage_cell(Particles,1);
Pinpart=zeros(size(particles,1),1);
flag=0;
for pci=1:1:numel(Particles)
    Pinpart(flag+1:flag+size(Particles{pci},1),1)=Part_list(pci);
    flag=flag+size(Particles{pci},1);
end
end