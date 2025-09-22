%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Nset,Sset,Eset]=Save_set_information(Set,Mesh)
Nset=Set.node_set;
Sset=Set.surface_set;
Eset=Set.element_set;
surfaces=Mesh.surfaces;
elements=Mesh.elements;
for nsi=1:1:numel(Sset)
    sid=Sset{nsi};
    Sset{nsi}=unique(surfaces(sid,:));
end
for nsi=1:1:numel(Eset)
    eid=Eset{nsi};
    nid=unique(elements(eid,1:end-1));
    Eset{nsi}=nid;
end
end