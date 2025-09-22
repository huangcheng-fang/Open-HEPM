%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Set=Resume_set_information(Nset,Sset,Eset,Mesh)
Anid=Mesh.new_nid;
node_num0=size(Mesh.nodes,1)-size(Anid,1);
for nsi=1:1:numel(Nset)
    nid=Nset{nsi};
    loc=ismember(Anid(:,1),nid)&ismember(Anid(:,end),nid);
    Nset{nsi}=[nid;find(loc)+node_num0];
end
surfaces=Mesh.surfaces;
for nsi=1:1:numel(Sset)
    nid=Sset{nsi};
    loc1=ismember(Anid(:,1),nid)&ismember(Anid(:,end),nid);
    nid0=[nid;find(loc1)+node_num0];
    loc2=ismember(surfaces(:,1),nid0)&ismember(surfaces(:,2),nid0);
    Sset{nsi}=find(loc2);
end
elements=Mesh.elements;
for nsi=1:1:numel(Eset)
    nid=Eset{nsi};
    loc1=ismember(Anid(:,1),nid)&ismember(Anid(:,end),nid);
    nid0=[nid;find(loc1)+node_num0];
    loc2=sum(ismember(elements(:,1:end-1),nid0),2)==size(elements,2)-1;
    Eset{nsi}=find(loc2);
end
Set.node_set=Nset;
Set.surface_set=Sset;
Set.element_set=Eset;
end